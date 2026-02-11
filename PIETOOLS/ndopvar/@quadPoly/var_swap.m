function P = var_swap(P, var1, var2)
%VAR_SWAP Swap one left variable's exponent with one right variable's exponent
% in a quadPoly stored with tensor-product monomial bases (cell Zs/Zt).
%
% If var1 is not in P.ns, it is appended as a constant-only variable (degree 0).
% If var2 is not in P.nt, it is appended as a constant-only variable (degree 0).
%
% After swapping, redundant degrees are merged and C is remapped accordingly.

    ns = P.ns; nt = P.nt;
    Zs = P.Zs; Zt = P.Zt;
    C  = sparse(P.C);

    m = P.dim(1);
    n = P.dim(2);

    % --- normalize empties ---
    if isempty(Zs), Zs = {}; end
    if isempty(Zt), Zt = {}; end
    if isempty(ns), ns = cell(1,numel(Zs)); end
    if isempty(nt), nt = cell(1,numel(Zt)); end

    % basis sizes
    ds = basisSize(Zs);
    dt = basisSize(Zt);

    % ---- ensure var1 exists in ns / Zs, var2 exists in nt / Zt ----
    idxS = find(strcmp(ns, var1), 1);
    if isempty(idxS)
        ns{end+1} = var1;
        Zs{end+1} = 0;           % constant-only degrees
        idxS = numel(ns);
        % update ds
        ds = basisSize(Zs);
    end

    idxT = find(strcmp(nt, var2), 1);
    if isempty(idxT)
        nt{end+1} = var2;
        Zt{end+1} = 0;           % constant-only degrees
        idxT = numel(nt);
        % update dt
        dt = basisSize(Zt);
    end

    % refresh locals after possible append
    P.ns = ns; P.nt = nt;
    P.Zs = Zs; P.Zt = Zt;

    % ---- if C is empty/zero, just commit the variable additions and return ----
    if isempty(C) || nnz(C)==0
        P.C = sparse(m*ds, n*dt);
        return;
    end

    % ---- find all nonzeros ----
    [r_idx, c_idx, val] = find(C);

    % Decode matrix row/col and full monomial indices
    rIn  = mod(r_idx - 1, m) + 1;     % 1..m
    cIn  = mod(c_idx - 1, n) + 1;     % 1..n
    kOld = ceil(r_idx / m);           % 1..ds (old full s-monomial)
    lOld = ceil(c_idx / n);           % 1..dt (old full t-monomial)

    % Unique (k,l) pairs
    pairKey = kOld + ds*(lOld-1);
    [uKey, ~, loc] = unique(pairKey);
    ku = mod(uKey - 1, ds) + 1;
    lu = floor((uKey - 1) / ds) + 1;

    % ---- convert full indices to per-variable position tuples ----
    sSizes = cellfun(@numel, Zs);
    tSizes = cellfun(@numel, Zt);

    sSubs = tensorDecomp(ku, sSizes);   % (#pairs)×|ns|
    tSubs = tensorDecomp(lu, tSizes);   % (#pairs)×|nt|

    % ---- convert positions -> actual degrees for each variable ----
    sDeg = subsToDegrees(sSubs, Zs);    % (#pairs)×|ns|
    tDeg = subsToDegrees(tSubs, Zt);    % (#pairs)×|nt|

    % ---- swap degrees between var1(left) and var2(right) ----
    tmp = sDeg(:, idxS);
    sDeg(:, idxS) = tDeg(:, idxT);
    tDeg(:, idxT) = tmp;

    % ---- rebuild per-variable degree lists to include only needed degrees ----
    [Zs_new, sPosMap] = rebuildBasisFromDegrees(Zs, sDeg);
    [Zt_new, tPosMap] = rebuildBasisFromDegrees(Zt, tDeg);

    ds_new = basisSize(Zs_new);
    dt_new = basisSize(Zt_new);

    % ---- map each unique pair to new full indices ----
    % old pair u -> new per-var positions from sPosMap/tPosMap applied to degrees
    sSubs_new = degreesToSubs(sDeg, Zs_new);
    tSubs_new = degreesToSubs(tDeg, Zt_new);

    kNew_u = subv2lin_rows(sSubs_new, cellfun(@numel, Zs_new)); % (#pairs)×1 in 1..ds_new
    lNew_u = subv2lin_rows(tSubs_new, cellfun(@numel, Zt_new)); % (#pairs)×1 in 1..dt_new

    % Map each nonzero entry via its pair id
    kNew_nz = kNew_u(loc);
    lNew_nz = lNew_u(loc);

    % ---- rebuild C (sparse sums duplicates automatically) ----
    r_new = (kNew_nz - 1) * m + rIn;
    c_new = (lNew_nz - 1) * n + cIn;

    C_new = sparse(r_new, c_new, val, m*ds_new, n*dt_new);

    % ---- commit ----
    P.ns = ns;
    P.nt = nt;
    P.Zs = Zs_new;
    P.Zt = Zt_new;
    P.C  = C_new;
end

% ====================== helpers ======================

function d = basisSize(Zcell)
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for i = 1:numel(Zcell)
        zi = Zcell{i};
        if isempty(zi), zi = 0; end
        d = d * numel(zi);
    end
end

function subs = tensorDecomp(idx, sizes)
% idx: vector in 1..prod(sizes); returns length(idx)×k subs (1-based)
    idx = idx(:) - 1;
    k = numel(sizes);
    subs = zeros(numel(idx), k);
    for i = k:-1:1
        si = sizes(i);
        subs(:,i) = mod(idx, si) + 1;
        idx = floor(idx / si);
    end
end

function deg = subsToDegrees(subs, Zcell)
% subs: L×k positions; returns L×k degrees
    L = size(subs,1);
    k = size(subs,2);
    deg = zeros(L,k);
    for i = 1:k
        zi = Zcell{i}(:);
        if isempty(zi), zi = 0; end
        deg(:,i) = zi(subs(:,i));
    end
end

function [Znew, posMap] = rebuildBasisFromDegrees(Zold, degUsed)
% For each variable i, set Znew{i} = unique(degUsed(:,i),'sorted')
% posMap{i}: map from degree value -> position in Znew{i} (via containers.Map avoided)
    k = numel(Zold);
    Znew = cell(1,k);
    posMap = cell(1,k); %#ok<NASGU>
    for i = 1:k
        z = unique(degUsed(:,i), 'sorted');
        if isempty(z), z = 0; end
        Znew{i} = z;
    end
end

function subs = degreesToSubs(deg, Zcell)
% deg: L×k degrees -> subs: L×k positions in Zcell lists
    L = size(deg,1);
    k = size(deg,2);
    subs = zeros(L,k);
    for i = 1:k
        zi = Zcell{i}(:);
        % degrees were built from zi so ismember should succeed
        [tf, loc] = ismember(deg(:,i), zi);
        if any(~tf)
            error('var_swap:degreeMissing', 'A needed degree is missing in rebuilt basis.');
        end
        subs(:,i) = loc;
    end
end

function idx = subv2lin_rows(subs, sizes)
% subs: L×k (1-based). sizes: 1×k. return L×1 linear indices (1..prod(sizes))
    L = size(subs,1);
    k = numel(sizes);
    idx0 = zeros(L,1);
    stride = ones(1,k);
    for i = 1:k-1
        stride(i) = prod(sizes(i+1:end));
    end
    stride(k) = 1;

    for i = 1:k
        idx0 = idx0 + (subs(:,i)-1) * stride(i);
    end
    idx = idx0 + 1;
end
