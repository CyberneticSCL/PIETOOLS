function obj = clean(obj, tol)
%clean Remove tiny coefficients, unused degrees, and constant-only variables.
%
%   obj = obj.clean(tol) or obj = clean(obj,tol) 
% Inputs:
%   obj: quadPoly object
%   tol (optional):
%       absolute threshold on |C_ij| below which entries are dropped.
%       If omitted, uses 1e-12 * max(abs(C(:))) (or 0 if C is all zeros).
%
% This reduction is *structure-preserving*: it only shrinks per-variable
% degree lists Zs{i}, Zt{j}. It does not remove arbitrary full monomials,
% because that would destroy the tensor-product basis structure.

    if nargin < 2 
        if isempty(obj.C) || nnz(obj.C)==0
            tol = 0;
        else
            tol = 1e-12*full(max(abs(obj.C(:))));
        end
    end

    m = obj.dim(1);
    n = obj.dim(2);

    % ---- 1) Drop small entries in C ----
    if ~isempty(obj.C) && nnz(obj.C) > 0 && tol > 0
        obj.C(abs(obj.C) < tol) = 0;
        obj.C = sparse(obj.C); % keep sparse tidy
    end

    % Early exit if everything is zero
    if isempty(obj.C) || nnz(obj.C) == 0
        % Minimal canonical "zero" object: keep 0-degree bases only
        obj.Zs = {}; obj.Zt = {};
        obj.ns = {}; obj.nt = {};
        obj.C = sparse(m, n);
        return;
    end

    % ---- 2) Determine which full s/t monomials are used (blockwise) ----
    ds_old = prod(cellfun(@numel, obj.Zs));
    dt_old = prod(cellfun(@numel, obj.Zt));

    % Used full monomials in s: any nonzero in corresponding m-row block
    usedS_full = false(ds_old,1);
    usedT_full = false(dt_old,1);

    % Find nonzeros once
    [rIdx, cIdx, ~] = find(obj.C);

    sMon = ceil(rIdx / m);
    tMon = ceil(cIdx / n);

    usedS_full(sMon) = true;
    usedT_full(tMon) = true;

    % ---- 3) Shrink per-variable degree sets, structure-preserving ----
    [Zs_new, ns_new, mapS_full] = reduceOneSide(obj.Zs, obj.ns, usedS_full);
    [Zt_new, nt_new, mapT_full] = reduceOneSide(obj.Zt, obj.nt, usedT_full);

    ds_new = prod(cellfun(@numel,Zs_new));
    dt_new = prod(cellfun(@numel,Zt_new));

    % ---- 4) Remap sparse coefficient matrix into reduced bases ----
    % Each entry C(row,col) belongs to:
    %   s full monomial = ceil(row/m), row inside block = mod(row-1,m)+1
    %   t full monomial = ceil(col/n), col inside block = mod(col-1,n)+1

    rIn = mod(rIdx-1, m) + 1;
    cIn = mod(cIdx-1, n) + 1;

    sMon_old = ceil(rIdx / m);
    tMon_old = ceil(cIdx / n);

    sMon_new = mapS_full(sMon_old);
    tMon_new = mapT_full(tMon_old);

    new_rIdx = (sMon_new - 1) * m + rIn;
    new_cIdx = (tMon_new - 1) * n + cIn;

    vals = full(obj.C(sub2ind(size(obj.C), rIdx, cIdx)));
    obj.C = sparse(new_rIdx, new_cIdx, vals, m*ds_new, n*dt_new);

    % ---- 5) Assign reduced bases/names and remove constant-only variables ----
    [obj.Zs, obj.ns] = dropConstantVars(Zs_new, ns_new);
    [obj.Zt, obj.nt] = dropConstantVars(Zt_new, nt_new);
end

% ========================= Local helpers =========================
function [Znew, nnew, mapFull] = reduceOneSide(Zcell, names, usedFull)
% Reduce per-variable degree lists based on which full monomials are used.
% Returns:
%   Znew: reduced Zcell
%   nnew: reduced names
%   mapFull: mapping from old full index -> new full index (0 if dropped)

    if isempty(Zcell)
        % No variables: only one monomial
        Znew = {};
        nnew = {};
        mapFull = 1; %#ok<NASGU>
        mapFull = ones(numel(usedFull),1);
        return;
    end

    k = numel(Zcell);
    dims_old = zeros(1,k);
    for i = 1:k, dims_old(i) = numel(Zcell{i}); end

    % Determine which positions in each variable are used
    usedPos = cell(1,k);
    for i = 1:k
        usedPos{i} = false(dims_old(i),1);
    end

    idxList = find(usedFull);
    for q = 1:numel(idxList)
        idx = idxList(q);
        subs = lin2subv(idx, dims_old);
        for i = 1:k
            usedPos{i}(subs(i)) = true;
        end
    end

    % Build reduced Z cell
    Znew = cell(1,k);
    nnew = names;
    posMap = cell(1,k); % old position -> new position
    dims_new = zeros(1,k);

    for i = 1:k
        keep = find(usedPos{i});
        if isempty(keep)
            % if nothing used, keep degree 0 as a safe constant
            keep = 1;
        end
        Znew{i} = Zcell{i}(keep);
        dims_new(i) = numel(keep);

        pm = zeros(dims_old(i),1);
        pm(keep) = (1:numel(keep)).';
        posMap{i} = pm;
    end

    % Map old full indices -> new full indices (0 if dropped)
    mapFull = zeros(numel(usedFull),1);
    idxList = find(usedFull);
    for q = 1:numel(idxList)
        idx = idxList(q);
        subs_old = lin2subv(idx, dims_old);

        subs_new = zeros(size(subs_old));
        for i = 1:k
            subs_new(i) = posMap{i}(subs_old(i));
            if subs_new(i) == 0
                % should not happen because we built keep from usedFull
                subs_new(i) = 0;
            end
        end
        mapFull(idx) = subv2lin(subs_new, dims_new);
    end
end

function subs = lin2subv(idx, dims)
% Inverse of subv2lin, with last dimension varying fastest.
    k = numel(dims);
    subs = zeros(1,k);
    r = idx - 1;
    for i = 1:k-1
        stride = prod(dims(i+1:end));
        subs(i) = floor(r / stride) + 1;
        r = r - (subs(i)-1)*stride;
    end
    subs(k) = r + 1;
end

function idx = subv2lin(subs, dims)
% Linear index from subscript vector, last dimension varies fastest.
    k = numel(dims);
    idx = 1;
    for i = 1:k-1
        stride = prod(dims(i+1:end));
        idx = idx + (subs(i)-1)*stride;
    end
    idx = idx + (subs(k)-1);
end

function [Zout, nout] = dropConstantVars(Zin, nin)
% Drop variables whose degree list is exactly [0] (constant-only).
    if isempty(Zin)
        Zout = {};
        nout = {};
        return;
    end
    keep = true(1, numel(Zin));
    for i = 1:numel(Zin)
        zi = Zin{i}(:);
        if zi==0
            keep(i) = false;
        end
    end
    Zout = Zin(keep);
    if isempty(nin)
        nout = cell(1, numel(Zout));
    else
        nout = nin(keep);
    end
end
