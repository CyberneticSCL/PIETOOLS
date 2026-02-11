function obj = combine(obj)
%COMBINE  Canonicalize a quadPoly by:
%  (1) sorting variable names (ns, nt) and reordering Zs/Zt accordingly,
%  (2) removing duplicate degrees within each per-variable exponent list,
%      and summing the corresponding coefficients in C.
%
% This preserves the quadPoly tensor-product basis structure:
%   Z(s) = kron_i s^{Zs{i}},   Z(t) = kron_j t^{Zt{j}}.

    % ---------- normalize empties ----------
    if isempty(obj.Zs), obj.Zs = {}; end
    if isempty(obj.Zt), obj.Zt = {}; end
    if isempty(obj.ns), obj.ns = cell(1, numel(obj.Zs)); end
    if isempty(obj.nt), obj.nt = cell(1, numel(obj.Zt)); end

    m = obj.dim(1);
    n = obj.dim(2);

    % ---------- build maps for left (s) and right (t) sides ----------
    [Zs_new, ns_new, mapS, ds_new] = combine_one_side(obj.Zs, obj.ns);
    [Zt_new, nt_new, mapT, dt_new] = combine_one_side(obj.Zt, obj.nt);

    % ---------- remap C using maps and sum duplicates ----------
    if isempty(obj.C) || nnz(obj.C) == 0
        obj.Zs = Zs_new; obj.ns = ns_new;
        obj.Zt = Zt_new; obj.nt = nt_new;
        obj.C  = sparse(m*ds_new, n*dt_new);
        return;
    end

    C = sparse(obj.C);
    [r, c, v] = find(C);

    % row block (matrix row), and monomial index in s-basis
    rIn    = mod(r-1, m) + 1;
    sMon   = ceil(r / m);
    sMon2  = mapS(sMon);

    % col block (matrix col), and monomial index in t-basis
    cIn    = mod(c-1, n) + 1;
    tMon   = ceil(c / n);
    tMon2  = mapT(tMon);

    % New global row/col indices in coefficient matrix
    r2 = (sMon2 - 1) * m + rIn;
    c2 = (tMon2 - 1) * n + cIn;

    % Sparse constructor sums duplicates automatically
    C2 = sparse(r2, c2, v, m*ds_new, n*dt_new);

    % ---------- write back ----------
    obj.Zs = Zs_new; obj.ns = ns_new;
    obj.Zt = Zt_new; obj.nt = nt_new;
    obj.C  = C2;

end

% ======================================================================
% Helper: combine one side (either s or t)
% ======================================================================
function [Z_new, names_new, mapFull, d_new] = combine_one_side(Z, names)
% Returns:
%   Z_new, names_new : sorted-by-name variables, degrees unique/sorted
%   mapFull          : mapping old full tensor index -> new full tensor index
%   d_new            : new basis size

    % No variables => basis size 1
    if isempty(Z)
        Z_new = {};
        names_new = {};
        mapFull = 1;
        d_new = 1;
        return;
    end

    k = numel(Z);

    % Ensure names exist
    if isempty(names)
        names = cell(1,k);
        for i = 1:k, names{i} = ''; end
    end

    % --------- (A) Sort variables by name (stable for ties) ----------
    % Keep unnamed variables at the end (but stable among themselves)
    key = names;
    for i = 1:k
        if isempty(key{i})
            key{i} = char(127); % high ASCII so empties go last
        end
    end
    [~, perm] = sort(key);
    Zsrt = Z(perm);
    nsrt = names(perm);

    sizes_srt = cellfun(@numel, Zsrt);
    d_srt = prod(sizes_srt);

    % Map old full index -> sorted-order full index
    % Step: old subs (in old order) -> subs in sorted order -> new linear index
    sizes_old = cellfun(@numel, Z);

    mapPerm = zeros(prod(sizes_old), 1);
    if isempty(mapPerm)  % if any size is 0 (shouldn't happen), treat as constant-only
        Z_new = {};
        names_new = {};
        mapFull = 1;
        d_new = 1;
        return;
    end

    for idx = 1:prod(sizes_old)
        subs_old = lin2subv(idx, sizes_old);      % 1×k, old order
        subs_srt = subs_old(perm);                % reorder to sorted variable order
        mapPerm(idx) = subv2lin(subs_srt, sizes_srt);
    end

    % --------- (B) Unique degrees within each variable (in sorted order) ----------
    Z_new = cell(1,k);
    names_new = nsrt;

    posMap = cell(1,k);  % posMap{i}(oldPos) = newPos for that variable
    sizes_new = zeros(1,k);

    for i = 1:k
        zi = Zsrt{i}(:);
        if isempty(zi)
            zi = 0;
        end
        [zu, ~, ic] = unique(zi, 'sorted');   % ic maps old positions -> new positions
        Z_new{i} = zu;
        sizes_new(i) = numel(zu);

        pm = ic(:); % 1..numel(zu)
        posMap{i} = pm;
    end

    d_new = prod(sizes_new);

    % Map sorted full index -> reduced full index
    mapDedup = zeros(d_srt, 1);
    for idx = 1:d_srt
        subs = lin2subv(idx, sizes_srt); % subs in sorted order
        subs2 = zeros(1,k);
        for i = 1:k
            subs2(i) = posMap{i}(subs(i));
        end
        mapDedup(idx) = subv2lin(subs2, sizes_new);
    end

    % --------- (C) Compose maps: old -> sorted -> dedup -----------
    mapFull = mapDedup(mapPerm);

    % --------- (D) Drop constant-only variables (optional, but usually desired) ----------
    % If a variable’s degree list is exactly [0], it contributes nothing to basis size.
    % Removing it is safe and makes a cleaner canonical form.
    keep = true(1,k);
    for i = 1:k
        zi = Z_new{i}(:);
        if numel(zi) == 1 && zi(1) == 0
            keep(i) = false;
        end
    end
    Z_new = Z_new(keep);
    names_new = names_new(keep);

    % NOTE: Dropping constant-only variables does NOT change d_new (multiplies by 1),
    % so mapFull and d_new remain valid.

end

% ======================================================================
% Index helpers (kron order, last dimension varies fastest)
% ======================================================================
function subs = lin2subv(idx, dims)
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
    k = numel(dims);
    idx = 1;
    for i = 1:k-1
        stride = prod(dims(i+1:end));
        idx = idx + (subs(i)-1)*stride;
    end
    idx = idx + (subs(k)-1);
end
