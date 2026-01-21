function P = reduce(P)
    m = P.dim(1);
    n = P.dim(2);

    Zs = double(P.Zs);
    Zt = double(P.Zt);
    C  = sparse(P.C);

    ds = size(Zs,1);
    dt = size(Zt,1);

    if ~isequal(size(C), [m*ds, n*dt])
        error('reduce: size(C) must be (m*size(Zs,1)) x (n*size(Zt,1)).');
    end

    % ---- merge duplicates in Zs (rows) ----
    [Zs_u, ~, gs] = unique(Zs, 'rows');   % gs maps old row -> new row
    ds_u = size(Zs_u,1);

    if ds_u ~= ds
        % Build row-aggregation matrix As: (ds_u x ds), As(gs(k),k)=1
        As = sparse(gs, 1:ds, 1, ds_u, ds);

        % Apply to each m-block of rows: equivalent to (I_m ⊗ As) * C
        A_big = kron(speye(m), As);
        C = A_big * C;

        Zs = Zs_u;
        ds = ds_u;
    end

    % ---- merge duplicates in Zt (rows) ----
    [Zt_u, ~, gt] = unique(Zt, 'rows');   % gt maps old row -> new row
    dt_u = size(Zt_u,1);

    if dt_u ~= dt
        % Column-aggregation matrix At: (dt_u x dt), At(gt(l),l)=1
        At = sparse(gt, 1:dt, 1, dt_u, dt);

        % Apply to each n-block of columns: equivalent to C * (I_n ⊗ At') 
        B_big = kron(speye(n), At');
        C = C * B_big;

        Zt = Zt_u;
        dt = dt_u;
    end

    % ---- commit ----
    P.Zs = Zs;
    P.Zt = Zt;
    P.C  = sparse(C);
end
