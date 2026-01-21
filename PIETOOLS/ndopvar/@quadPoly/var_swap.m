function P = var_swap(P, var1, var2)
    ns = P.ns; nt = P.nt;

    Zs = double(P.Zs);
    Zt = double(P.Zt);
    C  = sparse(P.C);

    m = P.dim(1);
    n = P.dim(2);

    ds = size(Zs,1);
    dt = size(Zt,1);

    % ---- Build final (standard) variable lists: ns gets var1, nt gets var2 ----
    idxS = find(strcmp(ns, var1), 1);
    idxT = find(strcmp(nt, var2), 1);

    ns_new = ns;
    if isempty(idxS)
        ns_new{end+1} = var1;
        idxS_new = numel(ns_new);
        addedS = true;
    else
        idxS_new = idxS;
        addedS = false;
    end

    nt_new = nt;
    if isempty(idxT)
        nt_new{end+1} = var2;
        idxT_new = numel(nt_new);
        addedT = true;
    else
        idxT_new = idxT;
        addedT = false;
    end

    % ---- Repack using only (k,l) pairs that appear in C ----
    [r_idx, c_idx, val] = find(C);

    % Decode basis indices
    k = mod(r_idx - 1, ds) + 1;          % 1..ds
    l = mod(c_idx - 1, dt) + 1;          % 1..dt
    i = floor((r_idx - 1) / ds) + 1;     % 1..m
    j = floor((c_idx - 1) / dt) + 1;     % 1..n

    % Unique (k,l) pairs
    pairKey = k + ds*(l-1);              % 1..ds*dt
    [uKey, ~, loc] = unique(pairKey);

    ku = mod(uKey - 1, ds) + 1;
    lu = floor((uKey - 1) / ds) + 1;

    % Exponents for each unique pair
    A = Zs(ku, :);   % numPairs x |ns|
    B = Zt(lu, :);   % numPairs x |nt|

    % Extend to ns_new/nt_new if we appended missing variables
    if addedS
        A = [A, zeros(size(A,1), 1)]; % var1 exponent slot starts at 0 on left
    end
    if addedT
        B = [B, zeros(size(B,1), 1)]; % var2 exponent slot starts at 0 on right
    end

    % Swap exponent "slots": var1 (left) <-> var2 (right)
    tmp = A(:, idxS_new);
    A(:, idxS_new) = B(:, idxT_new);
    B(:, idxT_new) = tmp;

    % Merge duplicates (build new bases)
    [Zs_new, ~, kmap_u] = unique(A, 'rows');
    [Zt_new, ~, lmap_u] = unique(B, 'rows');

    ds_new = size(Zs_new, 1);
    dt_new = size(Zt_new, 1);

    % Map each nonzero entry through its unique-pair id
    k_new_nz = kmap_u(loc);
    l_new_nz = lmap_u(loc);

    r_new = (i - 1) * ds_new + k_new_nz;
    c_new = (j - 1) * dt_new + l_new_nz;

    C_new = sparse(r_new, c_new, val, m*ds_new, n*dt_new);

    % ---- commit ----
    P.ns = ns_new;
    P.nt = nt_new;
    P.Zs = double(Zs_new);
    P.Zt = double(Zt_new);
    P.C  = C_new;
end
