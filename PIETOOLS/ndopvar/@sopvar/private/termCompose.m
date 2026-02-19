function [gammaLinIdx, CparamCells] = termCompose(Aqp, Bqp, alphaIdx, betaIdx, varsInB, varsOutA, varsInA, domMap, paramsSize, varsS3A)
%TERMCOMPOSE  Compose one term of A with one term of B for PI operators.
%
% This function computes the contribution of a single summand from A and a
% single summand from B in the composition C = A*B.
%
% INPUTS
%   Aqp, Bqp    : quadPoly terms from A.params{...} and B.params{...}
%   alphaIdx    : cell subscripts (per S3 variable) selecting A's PI mode
%   betaIdx     : cell subscripts (per S3 variable) selecting B's PI mode
%   varsInB     : B.vars_in   (C input variables)
%   varsOutA    : A.vars_out  (C output variables)
%   varsInA     : A.vars_in   (must match B.vars_out)
%   domMap      : OPTIONAL containers.Map baseVarName -> [a,b]
%                 If omitted/empty, defaults to [0,1] for all variables.
%   paramsSize  : OPTIONAL size of params cell array, e.g. size(A.params).
%                 If omitted, assumes 3 in every S3 dimension.
%
% OUTPUTS
%   gammaLinIdx : vector of linear indices into C.params
%   CparamCells : cell array of quadPoly terms to be added into C.params{gammaLinIdx(k)}
%
% Notes:
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

    if nargin < 8 || isempty(domMap)
        domMap = [];
    end

    % -------------------------------
    % 0) Index-dimension bookkeeping
    % -------------------------------
    nDim = numel(alphaIdx);               % number of subscript dimensions used by ind2sub
    paramsSize = paramsSize(:).';
    if numel(paramsSize) < nDim
        paramsSize(end+1:nDim) = 1;
    else
        paramsSize = paramsSize(1:nDim);
    end

    if nargin < 10 || isempty(varsS3A)
        % fallback, but you REALLY want to pass A.vars_S3
        varsS3A = intersect(varsInA, varsOutA, 'stable');
    end
    varsS3A = varsS3A(:).';
    nS3 = numel(varsS3A);

    % Map the k-th S3 variable to the correct dimension index in paramsSize/alphaIdx/betaIdx.
    % This is the key fix for the 1x3 vs 3x1 case.
    s3DimPos = mapS3Dims(paramsSize, nS3);

    % Pass-through variables of A (must be common to A input and output)
    S3pass = intersect(varsInA, varsOutA, 'stable');

    % --------------------------------------------
    % 1) Rename ALL intermediate variables to eta*
    %    Intermediate vars are exactly A's input vars (which equal B's output vars)
    % --------------------------------------------
    midVars = varsInA(:).';                 % 's...' names
    tMid    = cellfun(@s2t,   midVars, 'UniformOutput', false);
    etaMid  = cellfun(@s2eta, midVars, 'UniformOutput', false);

    % Ensure etaMid doesn't collide with existing names
    allNames = unique([Aqp.ns(:); Aqp.nt(:); Bqp.ns(:); Bqp.nt(:); tMid(:)]);
    etaMid   = makeUniqueNames(etaMid, allNames);

    % Build a quick lookup: etaOf('s3') -> 'eta3'
    etaOf = containers.Map(midVars, etaMid);

    % Rename A's right-side dummy vars t(..) -> eta(..)
    A2 = qpRename(Aqp, 't', tMid, etaMid);

    % Rename B's left-side output vars s(..) -> eta(..)
    B2 = qpRename(Bqp, 's', midVars, etaMid);

    % --------------------------------------------
    % 2) Multiply kernels in polynomial form
    % --------------------------------------------
    polyA = qpToPoly(A2);
    polyB = qpToPoly(B2);
    polyP = polyMtimes(polyA, polyB);

    % --------------------------------------------
    % 3) Fully integrate out A's S1 variables (input-only vars of A)
    %    IMPORTANT: integrate wrt the eta-version, not the s-version.
    % --------------------------------------------
    fullVars = setdiff(varsInA, varsOutA, 'stable');
    for i = 1:numel(fullVars)
        v = fullVars{i};
        [a,b] = getDom(domMap, v);
        polyP = polyIntegrate(polyP, etaOf(v), a, b);
    end

    % --------------------------------------------
    % 4) For each S3 variable of A, integrate out eta(s_k)
    %    with PI bounds determined by alphaIdx/betaIdx at the mapped dimension.
    %    Keep full gamma-subscript vectors (length nDim) per branch.
    % --------------------------------------------
    polys     = {polyP};
    gammaSubs = {ones(1, nDim)};   % each entry is a numeric row vector of subscripts

    for k = 1:nS3
        sVar = varsS3A{k};
        d    = s3DimPos(k);        % <-- the crucial mapping

        % If this S3-dimension is inactive in params (size 1), it contributes nothing.
        if paramsSize(d) == 1
            continue;
        end

        % If sVar is not actually pass-through in A, skip (should coincide with paramsSize(d)==1)
        if ~any(strcmp(S3pass, sVar))
            continue;
        end

        aOp = alphaIdx{d};   % 1/2/3
        bOp = betaIdx{d};    % 1/2/3

        etaVar = etaOf(sVar);  % e.g. 'eta1'
        thVar  = s2t(sVar);    % e.g. 't1'

        [aDom, bDom] = getDom(domMap, sVar);

        newPolys     = {};
        newGammaSubs = {};

        for r = 1:numel(polys)
            [plist, glist] = composeOneVar(polys{r}, sVar, etaVar, thVar, aDom, bDom, aOp, bOp);

            for q = 1:numel(plist)
                g = gammaSubs{r};
                g(d) = glist(q);   % set gamma along the correct PARAM DIMENSION

                newPolys{end+1,1}     = plist{q}; %#ok<AGROW>
                newGammaSubs{end+1,1} = g;        %#ok<AGROW>
            end
        end

        polys     = newPolys;
        gammaSubs = newGammaSubs;
    end

    % Sanity: eta-vars must be gone now
    for r = 1:numel(polys)
        etaLeft = polys{r}.vars(startsWith(polys{r}.vars,'eta'));
        if ~isempty(etaLeft)
            error('termCompose: eta variables survived elimination: %s', strjoin(etaLeft, ', '));
        end
    end

    % --------------------------------------------
    % 5) Build quadPoly outputs and gammaLinIdx
    % --------------------------------------------
    % Output vars (left side) are varsOutA
    leftVarsC = varsOutA(:).';

    % Input vars (right side) are dummy versions of varsInB
    % (this matches your s/t/eta convention)
    rightVarsC = cellfun(@s2t, varsInB(:).', 'UniformOutput', false);

    nRes = numel(polys);
    gammaLinIdx = zeros(nRes,1);
    CparamCells = cell(nRes,1);

    for r = 1:nRes
        gammaCell = num2cell(gammaSubs{r});
        gammaLinIdx(r) = sub2ind(paramsSize, gammaCell{:});

        preferredLeft = leftVarsC;
        CparamCells{r} = polyToQuadPoly(polys{r}, leftVarsC, rightVarsC, preferredLeft);
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
%                quadPoly renaming helper
% =====================================================================
function qpOut = qpRename(qpIn, side, oldNames, newNames)
% Rename variable names in quadPoly variable lists (ns or nt).
% side='s' -> ns (left vars), side='t' -> nt (right vars)

    qpOut = qpIn;

    if isempty(oldNames)
        return;
    end

    switch lower(side)
        case 's'
            v = qpOut.ns;
        case 't'
            v = qpOut.nt;
        otherwise
            error('side must be ''s'' or ''t''.');
    end

    for i = 1:numel(oldNames)
        idx = find(strcmp(v, oldNames{i}), 1);
        if ~isempty(idx)
            v{idx} = newNames{i};
        end
    end

    switch lower(side)
        case 's'
            qpOut.ns = v;
        case 't'
            qpOut.nt = v;
    end
end

function newNames = makeUniqueNames(baseNames, usedNames)
% Ensure each name in baseNames is unique w.r.t. usedNames and previously chosen names.
    newNames = baseNames(:);
    used = usedNames(:);

    for i = 1:numel(newNames)
        cand = newNames{i};
        suffix = 0;
        while any(strcmp(used, cand)) || any(strcmp(newNames(1:i-1), cand))
            suffix = suffix + 1;
            cand = sprintf('%s_%d', newNames{i}, suffix);
        end
        newNames{i} = cand;
    end
end

% =====================================================================
%                Domain lookup helper
% =====================================================================
function [a,b] = getDom(domMap, varName)
% domMap keyed by normal spatial vars: 's...'
% If varName is 't...' or 'eta...' map it back to 's...'

    if isempty(domMap)
        a = 0; b = 1; return;
    end

    base = varName;

    if startsWith(base,'t')
        base = ['s' base(2:end)];
    elseif startsWith(base,'eta')
        base = ['s' base(4:end)];
    end

    if isa(domMap,'containers.Map') && isKey(domMap, base)
        ab = domMap(base);
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

function key = expKey(e)
    key = sprintf('%d,', e);
end

function e = keyToExp(key)
    e = sscanf(key, '%d,').';
end

function poly = polyAddTerm(poly, e, M)
    key = expKey(e);
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
%
% Handles the MATLAB 1-var artifact: params can be 3x1 or 1x3 (both ndims=2).
    nDim = numel(paramsSize);

    if nS3 == 1
        d = find(paramsSize > 1, 1, 'first'); % the "3" dimension
        if isempty(d), d = 1; end
        dimPos = d;
        return;
    end

    % Typical case: dimensions correspond directly
    if nDim == nS3
        dimPos = 1:nS3;
        return;
    end

    % If extra dims exist, strip leading/trailing singleton artifacts until counts match.
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
