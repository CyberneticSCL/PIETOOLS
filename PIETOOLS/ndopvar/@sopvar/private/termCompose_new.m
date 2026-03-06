function [gammaLinIdx, CparamCells] = termCompose_new(...
                Aqp, Bqp, alphaIdx, betaIdx, varsInB, varsOutA, varsInA, domMap, paramsSize, varsS3A)
% Input quadPoly parameters have the form:
% A = Z(s4,s2a,s3a)' * CA[alpha] * Z(t2a,t2b,t3a,t3b) : 
% L2[s2a,s2b,s3a,s3b] -> L2[s4,s2a,s3a]

% B = Z(s2a,s2b,s3a,s3b)' * CB[beta] * Z(t1,t3a,t3b)
% L2[s1,s3a,s3b] -> L2[s2a,s2b,s3a,s3b]

% we need to find
% C = A*B
% C = Z(s4,s2a,s3a)'* CC[gamma] * Z(t1,t3b,t3a)
% where we compute CC[gamma] by doing
% C = Z(s4,s2a,s3a)'* CA[alpha] * INT(Z(eta2a,eta2b,eta3a,eta3b)
%                   *Z(eta2a,eta2b,eta3a,eta3b)') *CB[beta] * Z(t1,t3a,t3b)
% L2[s1,s3b,s3a] -> L2[s4,s2a,s3a]
% Summary of what happens:
%   1) Rename intermediate variables (A input == B output) to unique eta*
%   2) Multiply A and B kernels (as polynomials)
%   3) Integrate out A-only input variables (full integrals)
%   4) Eliminate each shared S3 eta-variable using PI rules (may branch)
%   5) Convert results back to quadPoly and return gamma indices

    % --------- Defaults / normalize shapes ----------
    if nargin < 8  || isempty(domMap),     domMap = []; end
    if nargin < 9  || isempty(paramsSize), paramsSize = [1 1]; end
    if nargin < 10 || isempty(varsS3A)
        varsS3A = intersect(varsInA, varsOutA, 'stable');
    end

    varsInB  = varsInB(:).';
    varsOutA = varsOutA(:).';
    varsInA  = varsInA(:).';
    varsS3A  = varsS3A(:).';

    % --------- 0) Index bookkeeping (alpha/beta/gamma live in paramsSize) ----------
    paramsSize = paramsSize(:).';
    nDim = max([numel(paramsSize), numel(alphaIdx), numel(betaIdx)]);
    paramsSize(end+1:nDim) = 1;  % pad with singleton dims if needed

    aSub = ones(1, nDim);
    bSub = ones(1, nDim);
    if ~isempty(alphaIdx), aSub(1:numel(alphaIdx)) = [alphaIdx{:}]; end
    if ~isempty(betaIdx),  bSub(1:numel(betaIdx))  = [betaIdx{:}];  end

    nS3      = numel(varsS3A);
    s3DimPos = mapS3Dims(paramsSize, nS3);  % which param-dim corresponds to each S3 var?
    S3pass   = intersect(varsInA, varsOutA, 'stable'); % pass-through vars of A

    % --------- 1) Rename intermediate vars (A input == B output) to eta* ----------
    midS   = varsInA;  % intermediate spatial vars (s...)
    midT   = cellfun(@s2t,   midS, 'UniformOutput', false);
    midEta = cellfun(@s2eta, midS, 'UniformOutput', false);

    if any(ismember(midEta, [Aqp.ns(:);Aqp.nt(:);Bqp.ns(:);Bqp.nt(:);midT(:)]))
        error('termCompose_new: eta collision.');
    end

    % A: rename t(midS) -> eta(midS)
    A2 = Aqp;  
    [tf,loc] = ismember(A2.nt, midT); 
    A2.nt(tf) = midEta(loc(tf));

    % B: rename s(midS) -> eta(midS)
    B2 = Bqp;  
    [tf,loc] = ismember(B2.ns, midS); 
    B2.ns(tf) = midEta(loc(tf));

    % Fast lookup: for any sVar in midS, etaVar is midEta{loc}
    [tfS3, locS3] = ismember(varsS3A, midS);
    if any(~tfS3)
        error('termCompose_new: varsS3A must be subset of varsInA.');
    end

    % --------- 2) Multiply kernels in polynomial form ----------
    poly = polyMtimes(qpToPoly(A2), qpToPoly(B2));

    % --------- 3) Fully integrate out A-only input variables ----------
    fullS = setdiff(varsInA, varsOutA, 'stable');  % input-only vars of A
    if ~isempty(fullS)
        [tfFull, locFull] = ismember(fullS, midS);
        if any(~tfFull)
            error('termCompose_new: fullS must be subset of varsInA.');
        end

        for i = 1:numel(fullS)
            sVar  = fullS{i};
            etaVar = midEta{locFull(i)};
            [a,b] = getDom(domMap, sVar);
            poly  = polyIntegrate(poly, etaVar, a, b);
        end
    end

    % --------- 4) Eliminate eta-vars for each shared S3 variable (branch if needed) ----------
    polyList  = {poly};
    gammaSubs = ones(1, nDim);   % one row per poly in polyList

    for k = 1:nS3
        sVar = varsS3A{k};

        d = s3DimPos(k);                 % param dimension for this S3 var
        if paramsSize(d) == 1, continue; end
        if ~any(strcmp(S3pass, sVar)), continue; end

        aOp = aSub(d);                   % 1/2/3
        bOp = bSub(d);                   % 1/2/3

        etaVar = midEta{locS3(k)};
        thVar  = s2t(sVar);
        [aDom, bDom] = getDom(domMap, sVar);

        nCur = numel(polyList);
        newPoly  = cell(2*nCur, 1);      % each branch can split into at most 2
        newGamma = zeros(2*nCur, nDim);

        pos = 0;
        for r = 1:nCur
            [plist, glist] = composeOneVar(polyList{r}, sVar, etaVar, thVar, aDom, bDom, aOp, bOp);

            nb = numel(plist);
            newPoly(pos+(1:nb)) = plist(:);

            gblock = repmat(gammaSubs(r,:), nb, 1);
            gblock(:, d) = glist(:);
            newGamma(pos+(1:nb), :) = gblock;

            pos = pos + nb;
        end

        polyList  = newPoly(1:pos);
        gammaSubs = newGamma(1:pos, :);
    end

    % Optional but helpful sanity check
    for r = 1:numel(polyList)
        if any(startsWith(polyList{r}.vars,'eta'))
            error('termCompose_new: eta variables survived elimination.');
        end
    end

    % --------- 5) Convert to quadPoly + compute linear gamma indices ----------
    leftVarsC  = varsOutA;
    rightVarsC = cellfun(@s2t, varsInB, 'UniformOutput', false);

    stride = [1, cumprod(paramsSize(1:end-1))];   % column-major strides
    gammaLinIdx = 1 + (gammaSubs - 1) * stride(:);

    nRes = numel(polyList);
    CparamCells = cell(nRes, 1);
    for r = 1:nRes
        CparamCells{r} = polyToQuadPoly(polyList{r}, leftVarsC, rightVarsC, leftVarsC);
    end
end

% =====================================================================
%                Composition for ONE shared variable
% =====================================================================
function [polyList, gammaList] = composeOneVar(polyIn, sVar, etaVar, thVar, aDom, bDom, aOp, bOp)
% Applies the 1D table of PI-mode combinations for variable (s, eta, theta).
%
% aOp, bOp, gamma are encoded as:
%   1 = M (multiplier)
%   2 = L (lower integral)
%   3 = U (upper integral)
%
% The returned polys have etaVar eliminated.

    % Multiplier on A collapses eta = s, gamma = bOp
    if aOp == 1
        polyList  = {polySubstitute(polyIn, etaVar, sVar)};
        gammaList = bOp;
        return;
    end

    % Multiplier on B collapses eta = theta, gamma = aOp
    if bOp == 1
        polyList  = {polySubstitute(polyIn, etaVar, thVar)};
        gammaList = aOp;
        return;
    end

    % Both are integrals: decide bounds for eta and possible branching
    if aOp == 2 && bOp == 2
        % (L,L) -> gamma=L, eta in [theta, s]
        polyList  = {polyIntegrate(polyIn, etaVar, thVar, sVar)};
        gammaList = 2;

    elseif aOp == 3 && bOp == 3
        % (U,U) -> gamma=U, eta in [s, theta]
        polyList  = {polyIntegrate(polyIn, etaVar, sVar, thVar)};
        gammaList = 3;

    elseif aOp == 2 && bOp == 3
        % (L,U) -> split:
        %   gamma=L : eta in [aDom, theta]
        %   gamma=U : eta in [aDom, s]
        polyList  = { ...
            polyIntegrate(polyIn, etaVar, aDom, thVar), ...
            polyIntegrate(polyIn, etaVar, aDom, sVar) ...
        };
        gammaList = [2, 3];

    elseif aOp == 3 && bOp == 2
        % (U,L) -> split:
        %   gamma=L : eta in [s, bDom]
        %   gamma=U : eta in [theta, bDom]
        polyList  = { ...
            polyIntegrate(polyIn, etaVar, sVar,  bDom), ...
            polyIntegrate(polyIn, etaVar, thVar, bDom) ...
        };
        gammaList = [2, 3];

    else
        error('Unexpected PI-mode combination (aOp=%d, bOp=%d).', aOp, bOp);
    end
end

% =====================================================================
%                Domain lookup helper
% =====================================================================
function [a,b] = getDom(domMap, sVar)
% domMap keyed by normal spatial vars: 's...'
    if isempty(domMap)
        a = 0; b = 1; return;
    end
    if isa(domMap,'containers.Map') && isKey(domMap, sVar)
        ab = domMap(sVar);
        a = ab(1); b = ab(2);
    else
        a = 0; b = 1;
    end
end

% =====================================================================
%      Polynomial "dictionary" structure and operations
% =====================================================================

function poly = polyCreate(vars, m, n)
    poly.vars = vars(:)'; % row cell
    poly.m = m;
    poly.n = n;
    poly.map = containers.Map('KeyType','char','ValueType','any');
end


function e = keyToExp(key)
    e = sscanf(key, '%d,').';
end

function poly = polyAddTerm(poly, e, M)
    key = sprintf('%d,', e);
    if isKey(poly.map, key)
        poly.map(key) = poly.map(key) + M;
    else
        poly.map(key) = M;
    end
end

function polyOut = polyReorder(polyIn, newVars)
% Reorder/extend polynomial variables to match newVars.
    polyOut = polyCreate(newVars, polyIn.m, polyIn.n);

    [tf, loc] = ismember(polyIn.vars, newVars);
    if any(~tf)
        error('polyReorder: some variables are missing from newVars.');
    end

    keys = polyIn.map.keys;
    for i = 1:numel(keys)
        eOld = keyToExp(keys{i});
        M = polyIn.map(keys{i});

        eNew = zeros(1, numel(newVars));
        eNew(loc) = eOld;

        polyOut = polyAddTerm(polyOut, eNew, M);
    end
end

function polyC = polyMtimes(polyA, polyB)
% Polynomial-matrix product: sum over monomials with exponent addition and matrix multiply.

    if polyA.n ~= polyB.m
        error('polyMtimes: inner dimensions do not match.');
    end

    varsU = union(polyA.vars, polyB.vars, 'stable');
    A = polyReorder(polyA, varsU);
    B = polyReorder(polyB, varsU);

    polyC = polyCreate(varsU, polyA.m, polyB.n);

    keysA = A.map.keys;
    keysB = B.map.keys;

    for i = 1:numel(keysA)
        eA = keyToExp(keysA{i});
        MA = A.map(keysA{i});

        for j = 1:numel(keysB)
            eB = keyToExp(keysB{j});
            MB = B.map(keysB{j});

            eC = eA + eB;
            MC = MA * MB;

            polyC = polyAddTerm(polyC, eC, MC);
        end
    end
end

function polyOut = polySubstitute(polyIn, varX, bound)
% Substitute varX = bound (bound can be numeric or a variable name).
% This is used for multiplier (delta) collapses: no integration scaling.

    idxX = find(strcmp(polyIn.vars, varX), 1);
    if isempty(idxX)
        polyOut = polyIn;
        return;
    end

    varsNew = polyIn.vars;
    varsNew(idxX) = [];
    polyOut = polyCreate(varsNew, polyIn.m, polyIn.n);

    keys = polyIn.map.keys;
    for i = 1:numel(keys)
        e = keyToExp(keys{i});
        M = polyIn.map(keys{i});

        k = e(idxX);
        e(idxX) = [];
        if isnumeric(bound)
            Mnew = M * (bound^k);
            enew = e;
        else
            idxYold = find(strcmp(polyIn.vars, bound), 1);
            if isempty(idxYold)
                error('polySubstitute: bound variable %s not found.', bound);
            end
            idxY = idxYold - (idxYold > idxX);
            enew = e;
            enew(idxY) = enew(idxY) + k;
            Mnew = M;
        end

        polyOut = polyAddTerm(polyOut, enew, Mnew);
    end
end

function polyOut = polyIntegrate(polyIn, varX, lower, upper)
% Compute integral ∫_{lower}^{upper} polyIn d(varX)
% where lower/upper can be numeric or a variable name.
%
% This eliminates varX.

    idxX = find(strcmp(polyIn.vars, varX), 1);
    if isempty(idxX)
        error('polyIntegrate: variable %s not present (unexpected in this usage).', varX);
    end

    varsNew = polyIn.vars;
    varsNew(idxX) = [];
    polyOut = polyCreate(varsNew, polyIn.m, polyIn.n);

    keys = polyIn.map.keys;
    for i = 1:numel(keys)
        e0 = keyToExp(keys{i});
        M0 = polyIn.map(keys{i});

        k = e0(idxX);
        scale = 1/(k+1);

        % Remove varX from exponent vector
        e = e0;
        e(idxX) = [];

        % Upper contribution
        [eU, MU] = evalAntiderivativeAt(e0, M0, idxX, upper, +scale);
        polyOut = polyAddTerm(polyOut, eU, MU);

        % Lower contribution (subtract)
        [eL, ML] = evalAntiderivativeAt(e0, M0, idxX, lower, -scale);
        polyOut = polyAddTerm(polyOut, eL, ML);
    end

    function [eNew, MNew] = evalAntiderivativeAt(eFull, MFull, idxXfull, boundVal, signedScale)
        kLocal = eFull(idxXfull);

        % exponent vector without x
        eNoX = eFull;
        eNoX(idxXfull) = [];

        if isnumeric(boundVal)
            % eta^(k+1) -> bound^(k+1)
            MNew = MFull * (boundVal^(kLocal+1)) * signedScale;
            eNew = eNoX;
        else
            % eta^(k+1) -> (boundVar)^(k+1) i.e. shift exponent to boundVar
            idxYold = find(strcmp(polyIn.vars, boundVal), 1);
            if isempty(idxYold)
                error('polyIntegrate: bound variable %s not found.', boundVal);
            end
            idxY = idxYold - (idxYold > idxXfull);

            eNew = eNoX;
            eNew(idxY) = eNew(idxY) + (kLocal+1);

            MNew = MFull * signedScale;
        end
    end
end

% =====================================================================
%   quadPoly <-> polynomial dictionary conversion
% =====================================================================

function poly = qpToPoly(qp)
% Convert quadPoly to a monomial dictionary poly(vars -> coeff matrix).
%
% The result is:
%   poly(vars) = sum_{e} M_e * prod_i vars{i}^{e(i)}

    vars = [qp.ns(:); qp.nt(:)].';
    m = qp.dim(1); n = qp.dim(2);

    poly = polyCreate(vars, m, n);

    % Handle empty Z lists
    if isempty(qp.Zs), ds = 1; expS = zeros(1,0); else, expS = basisExponents(qp.Zs); ds = size(expS,1); end
    if isempty(qp.Zt), dt = 1; expT = zeros(1,0); else, expT = basisExponents(qp.Zt); dt = size(expT,1); end

    if isempty(qp.C)
        return;
    end

    [rr, cc, vv] = find(qp.C);

    for idx = 1:numel(vv)
        r = rr(idx); c = cc(idx); val = vv(idx);

        outRow = floor((r-1)/ds) + 1;
        iS     = r - (outRow-1)*ds;

        inCol  = floor((c-1)/dt) + 1;
        iT     = c - (inCol-1)*dt;

        e = [expS(iS,:), expT(iT,:)];

        M = zeros(m,n);
        M(outRow, inCol) = val;

        poly = polyAddTerm(poly, e, M);
    end
end

function expMat = basisExponents(Zcell)
% Return exponent combinations for the Kronecker basis Z(s1) ⊗ Z(s2) ⊗ ...
% Ordering assumes "last variable changes fastest" (Kronecker product convention).

    nv = numel(Zcell);
    if nv == 0
        expMat = zeros(1,0);
        return;
    end

    lens = cellfun(@numel, Zcell);
    d = prod(lens);

    expMat = zeros(d, nv);

    for idx = 1:d
        tmp = idx - 1;
        subs = zeros(1, nv);
        for i = nv:-1:1
            L = lens(i);
            subs(i) = mod(tmp, L) + 1;
            tmp = floor(tmp / L);
        end
        for i = 1:nv
            expMat(idx,i) = Zcell{i}(subs(i));
        end
    end
end

function qp = polyToQuadPoly(poly, leftVars, rightVars, preferredLeft)
% Convert polynomial dictionary back to quadPoly with specified left/right split.
%
% If poly contains variables not listed in leftVars/rightVars, this function
% will automatically assign them to a side:
%   - if variable is in preferredLeft => left
%   - else if variable endswith '_dum' => right
%   - else => right (safe default for kernel "input" vars)

    if nargin < 4
        preferredLeft = {};
    end
    preferredLeft = preferredLeft(:)';

    % Ensure row-cells
    leftVars  = leftVars(:)';
    rightVars = rightVars(:)';

    % ------------------------------------------------------------------
    % 0) Ensure every poly variable appears in leftVars or rightVars
    % ------------------------------------------------------------------
    allLR = union(leftVars, rightVars, 'stable');
    missing = setdiff(poly.vars, allLR, 'stable');

    if ~isempty(missing)
        % Auto-assign missing variables
        addLeft  = {};
        addRight = {};

        for i = 1:numel(missing)
            v = missing{i};
            if any(strcmp(preferredLeft, v)) || startsWith(v,'s')
                addLeft{end+1} = v;
            elseif startsWith(v,'t')
                addRight{end+1} = v;
            elseif startsWith(v,'eta')
                error('polyToQuadPoly: eta variable survived elimination: %s', v);
            else
                % safe default: put on right
                addRight{end+1} = v;
            end
        end

        leftVars  = [leftVars,  addLeft];
        rightVars = [rightVars, addRight];
    end

    % Remove any accidental duplicates while preserving order
    leftVars  = unique(leftVars, 'stable');
    rightVars = unique(rightVars, 'stable');

    % If something ended up on both sides, force it to the left (arbitrary but consistent)
    both = intersect(leftVars, rightVars, 'stable');
    if ~isempty(both)
        rightVars = setdiff(rightVars, both, 'stable');
    end

    % ------------------------------------------------------------------
    % 1) Reorder poly vars to [leftVars, rightVars]
    % ------------------------------------------------------------------
    allVars = [leftVars(:).', rightVars(:).'];
    poly = polyReorder(poly, allVars);

    nL = numel(leftVars);
    nR = numel(rightVars);

    % ------------------------------------------------------------------
    % 2) Collect exponent sets per variable
    % ------------------------------------------------------------------
    exps = cell(1, numel(allVars));
    for i = 1:numel(allVars)
        exps{i} = 0; % always include 0
    end

    keys = poly.map.keys;
    for i = 1:numel(keys)
        e = keyToExp(keys{i});
        for v = 1:numel(allVars)
            exps{v} = unique([exps{v}(:); e(v)]);
        end
    end
    for v = 1:numel(allVars)
        exps{v} = sort(exps{v}(:));
    end

    Zs = exps(1:nL);
    Zt = exps(nL+1:end);

    ds = prod(cellfun(@numel, Zs));
    dt = prod(cellfun(@numel, Zt));

    m = poly.m; n = poly.n;

    I = zeros(0,1); J = zeros(0,1); V = zeros(0,1);

    lenL = cellfun(@numel, Zs);
    lenR = cellfun(@numel, Zt);

    % ------------------------------------------------------------------
    % 3) Fill sparse C using Kronecker indices
    % ------------------------------------------------------------------
    for i = 1:numel(keys)
        e = keyToExp(keys{i});
        M = poly.map(keys{i});

        eL = e(1:nL);
        eR = e(nL+1:end);

        posL = zeros(1,nL);
        for k = 1:nL
            posL(k) = find(Zs{k} == eL(k), 1);
        end
        posR = zeros(1,nR);
        for k = 1:nR
            posR(k) = find(Zt{k} == eR(k), 1);
        end

        idxL = kronIndex(posL, lenL);
        idxR = kronIndex(posR, lenR);

        [ir, ic, val] = find(M);
        for t = 1:numel(val)
            row = (ir(t)-1)*ds + idxL;
            col = (ic(t)-1)*dt + idxR;
            I(end+1,1) = row; %#ok<AGROW>
            J(end+1,1) = col; %#ok<AGROW>
            V(end+1,1) = val(t); %#ok<AGROW>
        end
    end

    C = sparse(I, J, V, m*ds, n*dt);

    qp = quadPoly(C, Zs, Zt, [m,n], leftVars(:).', rightVars(:).');
end


function idx = kronIndex(pos, lens)
% Convert per-variable basis positions (1-based) to a single Kronecker index (1-based)
% with "last variable changes fastest".
    if isempty(lens)
        idx = 1;
        return;
    end
    idx = 1;
    stride = 1;
    for i = numel(lens):-1:1
        idx = idx + (pos(i)-1)*stride;
        stride = stride * lens(i);
    end
end
function t = s2t(s)
% Map 'sXYZ' -> 'tXYZ'
    if startsWith(s,'s')
        t = ['t' s(2:end)];
    else
        error('s2t expects a variable starting with ''s'': %s', s);
    end
end

function e = s2eta(s)
% Map 'sXYZ' -> 'etaXYZ'
    if startsWith(s,'s')
        e = ['eta' s(2:end)];
    else
        error('s2eta expects a variable starting with ''s'': %s', s);
    end
end

function dimPos = mapS3Dims(paramsSize, nS3)
% Map each S3 variable index to a dimension index in paramsSize.
% Handles MATLAB 1-var artifact: params can be 3x1 or 1x3.

    if nS3 == 0
        dimPos = [];
        return;
    end

    nDim = numel(paramsSize);

    if nS3 == 1
        d = find(paramsSize > 1, 1, 'first'); % the "3" dimension
        if isempty(d), d = 1; end
        dimPos = d;
        return;
    end

    if nDim == nS3
        dimPos = 1:nS3;
        return;
    end

    dims = 1:nDim;
    while numel(dims) > nS3 && paramsSize(dims(end)) == 1
        dims(end) = [];
    end
    while numel(dims) > nS3 && paramsSize(dims(1)) == 1
        dims(1) = [];
    end
    if numel(dims) ~= nS3
        error('mapS3Dims: cannot map %d S3 vars onto paramsSize=%s', nS3, mat2str(paramsSize));
    end
    dimPos = dims;
end
