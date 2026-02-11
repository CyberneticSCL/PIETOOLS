function Q = polynomial2quadPoly(P, sNames, tNames)
%POLYNOMIAL2QUADPOLY Convert polynomial object to quadPoly.
%
% Inputs:
%   P      : polynomial object
%   sNames : cellstr of left variables (ns)
%   tNames : cellstr of right variables (nt)
%
% Output:
%   Q      : quadPoly with tensor-product bases per variable.
%
% Notes:
%   - Builds Zs/Zt as unique degrees per variable seen in P.degmat.
%   - Requires that P.varname contains only sNames∪tNames (constants allowed).

    m = P.matdim(1);
    n = P.matdim(2);

    var = P.varname(:);
    deg = full(P.degmat);          % T×V (can be big but manageable)
    coeff = sparse(P.coefficient); % T×(m*n)

    % Map polynomial variables into left/right positions
    [tfS, posS] = ismember(sNames(:), var);
    [tfT, posT] = ismember(tNames(:), var);

    if any(~tfS)
        error('polynomial_to_quadPoly:missingVar', 'Some sNames not present in P.varname.');
    end
    if any(~tfT)
        error('polynomial_to_quadPoly:missingVar', 'Some tNames not present in P.varname.');
    end

    % Build Zs/Zt degree lists (unique degrees seen per variable)
    Zs = cell(1, numel(sNames));
    for i = 1:numel(sNames)
        Zs{i} = unique(deg(:, posS(i)), 'sorted');
        if isempty(Zs{i}), Zs{i} = 0; end
    end

    Zt = cell(1, numel(tNames));
    for i = 1:numel(tNames)
        Zt{i} = unique(deg(:, posT(i)), 'sorted');
        if isempty(Zt{i}), Zt{i} = 0; end
    end

    ds = basisSize(Zs);
    dt = basisSize(Zt);

    % Precompute maps from degree tuples -> tensor index for s and t
    % We use ismember per variable and then linearize subscripts.
    sSizes = cellfun(@numel, Zs);
    tSizes = cellfun(@numel, Zt);
    sStride = kronStrides(sSizes);
    tStride = kronStrides(tSizes);

    % For each term row in P, compute ks (1..ds) and kt (1..dt)
    T = size(deg,1);
    ks = ones(T,1);
    kt = ones(T,1);

    for i = 1:numel(sNames)
        [ok, loc] = ismember(deg(:, posS(i)), Zs{i});
        if any(~ok), error('polynomial_to_quadPoly:degreeMissingS','Degree missing in Zs at %s', sNames{i}); end
        ks = ks + (loc-1) * sStride(i);
    end

    for i = 1:numel(tNames)
        [ok, loc] = ismember(deg(:, posT(i)), Zt{i});
        if any(~ok), error('polynomial_to_quadPoly:degreeMissingT','Degree missing in Zt at %s', tNames{i}); end
        kt = kt + (loc-1) * tStride(i);
    end

    % Now place coefficients into Q.C:
    % Q.C row = (i-1)*ds + ks(term)
    % Q.C col = (j-1)*dt + kt(term)
    %
    % We build by accumulating sparse triplets.
    [termIdx, entryIdx, val] = find(coeff); % coeff(term, entry)

    % entryIdx corresponds to matrix entry col in polynomial format:
    % entry = i + m*(j-1)  => i = mod(entry-1,m)+1, j = floor((entry-1)/m)+1
    iMat = mod(entryIdx-1, m) + 1;
    jMat = floor((entryIdx-1)/m) + 1;

    row = (iMat-1) * ds + ks(termIdx);
    col = (jMat-1) * dt + kt(termIdx);

    C = sparse(row, col, val, m*ds, n*dt);

    % Construct quadPoly (property assignment is fine if no custom ctor)
    Q = quadPoly();
    Q.C = C;
    Q.Zs = Zs;
    Q.Zt = Zt;
    Q.dim = [m n];
    Q.ns = sNames(:).';
    Q.nt = tNames(:).';
end

% ---------------- helpers ----------------
function d = basisSize(Zcell)
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for k = 1:numel(Zcell)
        z = Zcell{k};
        if isempty(z), z = 0; end
        d = d * numel(z);
    end
end

function stride = kronStrides(sizes)
% last variable varies fastest
    k = numel(sizes);
    stride = ones(1,k);
    for i = 1:k-1
        stride(i) = prod(sizes(i+1:end));
    end
    stride(k) = 1;
end
