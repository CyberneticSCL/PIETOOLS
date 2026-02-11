function [Minv_qp, RM1_qp, RM2_qp, info] = ...
   inv_1D(M_qp, FG1_qp, FG2_qp, interval, opts)
% inv_1D
% Given [a,b], M, and FG1/FG2
% representing the operator
%   (M - K)x = M(s)x(s) - ∫_a^s M(s)^{-1} FG1(s,t) x(t) dt - ∫_s^b M(s)^{-1} FG2(s,t)x(t) dt
%
% Compute (M - K)^(-1) = (I - K)^(-1) M^{-1} = M^{-1} + R M^{-1}
% using Gohberg–Krein, where
%   K1(s,t) = M(s)^{-1} FG1(s,t) ≈ M(s)^{-1}F1(s)G1(t)
%   K2(s,t) = M(s)^{-1} FG2(s,t) ≈ M(s)^{-1}F2(s)G2(t)
%
% Inputs
%   M_qp   : quadPoly, m×m, depends only on s
%   FG1_qp : quadPoly kernel, m×m, is NOT split
%   FG2_qp : quadPoly kernel, m×m, is NOT split
%   interval : [a,b]
%   opts   : optional struct with fields:
%       .tol_rank (default 1e-12)  numerical rank tol for decomposeLR
%       .alpha    (default 3)      oversampling factor for basis -> grid size
%       .Nmin     (default 8)
%       .maxN     (default 512)
%       .rk4_cfl  (default 0.5)    enforce h*L <= rk4_cfl
%       .storeABC (default true)   store A(t) on grid for ODE (recommended)
%
% Outputs
%   Minv_qp : quadPoly approx of M(s)^(-1) (left-only)
%   RM1_qp  : quadPoly kernel for (R M^{-1})(s,t) for t <= s
%   RM2_qp  : quadPoly kernel for (R M^{-1})(s,t) for t >  s
%   info    : diagnostics (ranks, grid, rcond, etc.)
%
% Notes
%   - Lambda is NOT used inside. FG1_qp/FG2_qp must already include it.
%   - We solve U' = A(s) U, V' = -V A(s), where A(s)=B(s)C(s).
%   - Fits Minv and RM kernels back to quadPoly only at the end.

    % ----------------------- Defaults -----------------------
    if nargin < 5, opts = struct(); end
    if ~isfield(opts,'tol_rank'), opts.tol_rank = 1e-12; end
    if ~isfield(opts,'alpha'),    opts.alpha    = 3;     end
    if ~isfield(opts,'Nmin'),     opts.Nmin     = 8;     end
    if ~isfield(opts,'maxN'),     opts.maxN     = 512;   end
    if ~isfield(opts,'rk4_cfl'),  opts.rk4_cfl  = 0.5;   end
    if ~isfield(opts,'storeABC'), opts.storeABC = true;  end

    a = interval(1); b = interval(2);
    m = M_qp.dim(1);

    % ----------------- 0) Split FG1/FG2 first: FG ≈ F(s)G(t) -----------------
    [F1_qp, G1_qp, infoFG1] = decomposeLR(FG1_qp, opts.tol_rank);
    [F2_qp, G2_qp, infoFG2] = decomposeLR(FG2_qp, opts.tol_rank);

    n1 = F1_qp.dim(2);
    n2 = F2_qp.dim(2);
    n  = n1 + n2;

    if n == 0
        % K is numerically zero -> inverse is just M^{-1}
        % (compute_Minv_leftOnly not used here; keep it simple)
        dsM = quadPoly_basisSize(M_qp.Zs);
        N0  = min(opts.maxN, max(opts.Nmin, ceil(opts.alpha*dsM)));
        sgrid = linspace(a,b,N0);

        Minv = zeros(m,m,N0);
        for i = 1:N0
            Mi = quadPoly_eval(M_qp, sgrid(i), []);
            Minv(:,:,i) = Mi \ eye(m);
        end
        Minv_qp = fit_leftOnly_quadPoly(M_qp.Zs, M_qp.ns, sgrid, Minv);

        RM1_qp = zeroKernelLike(FG1_qp);
        RM2_qp = zeroKernelLike(FG1_qp);

        info = struct('note','FG ranks are zero; K≈0', 'N',N0, 'rank1',n1, 'rank2',n2);
        return;
    end

    % ----------------- 1) Choose grid size N -----------------
    ds_out = quadPoly_basisSize(FG1_qp.Zs);
    dt_in  = quadPoly_basisSize(FG1_qp.Zt);

    N_basis = max(opts.Nmin, ceil(opts.alpha * max(ds_out, dt_in)));
    N = min(opts.maxN, N_basis);
    sgrid = linspace(a, b, N);
    h = (b-a)/(N-1);

    % ----------------- 2) Evaluate everything needed on the grid -----------------
    Minv = zeros(m,m,N);
    C = zeros(m,n,N);
    B = zeros(n,m,N);
    if opts.storeABC
        A = zeros(n,n,N);
    else
        A = [];
    end

    for i = 1:N
        s = sgrid(i);

        Mi = quadPoly_eval(M_qp, s, []);
        Minv(:,:,i) = Mi \ eye(m);

        F1i = quadPoly_eval(F1_qp, s, []);          % m×n1
        F2i = quadPoly_eval(F2_qp, s, []);          % m×n2
        G1i = quadPoly_eval(G1_qp, [], s);          % n1×m
        G2i = quadPoly_eval(G2_qp, [], s);          % n2×m

        C(:,:,i) = [Minv(:,:,i)*F1i, Minv(:,:,i)*F2i]; % m×n
        B(:,:,i) = [G1i; -G2i];                         % n×m

        if opts.storeABC
            A(:,:,i) = B(:,:,i) * C(:,:,i);            % n×n
        end
    end

    % ----------------- 3) Optional CFL refinement using coarse A norm -----------------
    % NOTE: lambda removed, so heuristic is h*||A|| <= rk4_cfl
    Lest = 0;
    Nc = min(20, max(10, floor(N/4)));
    idxc = unique(round(linspace(1, N, Nc)));

    if opts.storeABC
        for k = 1:numel(idxc)
            Lest = max(Lest, norm(A(:,:,idxc(k)), 'fro'));
        end
    else
        for k = 1:numel(idxc)
            ii = idxc(k);
            Lest = max(Lest, norm(B(:,:,ii) * C(:,:,ii), 'fro'));
        end
    end

    N_cfl = 0;
    if Lest > 0
        N_cfl = ceil((b-a) * Lest / max(opts.rk4_cfl, eps));
    end

    N_new = min(opts.maxN, max([N, N_cfl, opts.Nmin]));
    if N_new ~= N
        opts2 = opts; opts2.Nmin = N_new; opts2.alpha = 1;
        [Minv_qp, RM1_qp, RM2_qp, info] = sopvar.inv_1D(M_qp, FG1_qp, FG2_qp, interval, opts2);
        return;
    end

    % ----------------- 4) Propagate U and V on the grid (RK4) -----------------
    % Now: U' = A U,   V' = -V A
    U = zeros(n,n,N);
    V = zeros(n,n,N);
    U(:,:,1) = eye(n);
    V(:,:,1) = eye(n);

    for i = 1:N-1
        Ai   = A_at(i);
        Aip1 = A_at(i+1);
        Ahalf = 0.5*(Ai + Aip1);

        % U' = A U
        k1U = Ai   * U(:,:,i);
        k2U = Ahalf * (U(:,:,i) + (h/2)*k1U);
        k3U = Ahalf * (U(:,:,i) + (h/2)*k2U);
        k4U = Aip1  * (U(:,:,i) + h*k3U);
        U(:,:,i+1) = U(:,:,i) + (h/6)*(k1U + 2*k2U + 2*k3U + k4U);

        % V' = -V A
        k1V = -V(:,:,i)              * Ai;
        k2V = -(V(:,:,i)+(h/2)*k1V)  * Ahalf;
        k3V = -(V(:,:,i)+(h/2)*k2V)  * Ahalf;
        k4V = -(V(:,:,i)+h*k3V)      * Aip1;
        V(:,:,i+1) = V(:,:,i) + (h/6)*(k1V + 2*k2V + 2*k3V + k4V);
    end

    % ----------------- 5) Build P from U(b) using a solve -----------------
    Ub = U(:,:,end);
    U21b = Ub(n1+1:end, 1:n1);
    U22b = Ub(n1+1:end, n1+1:end);

    info.rcondU22b = rcond(U22b);
    if info.rcondU22b < 1e-14
        warning('U22(b) appears near singular (rcond=%g). Invertibility may fail.', info.rcondU22b);
    end

    X = U22b \ U21b; % solve, not explicit inverse

    P = zeros(n,n);
    if n2 > 0
        P(n1+1:end, 1:n1) = X;
        P(n1+1:end, n1+1:end) = eye(n2);
    end
    IminusP = eye(n) - P;

    % ----------------- 6) Precompute T_i and S_j -----------------
    T1 = zeros(m,n,N);
    T2 = zeros(m,n,N);
    S  = zeros(n,m,N);

    for i = 1:N
        T1(:,:,i) = C(:,:,i) * U(:,:,i) * IminusP;
        T2(:,:,i) = -C(:,:,i) * U(:,:,i) * P;
        S(:,:,i)  = V(:,:,i) * B(:,:,i) * Minv(:,:,i);
    end

    % ----------------- 7) Fit Minv(s) to quadPoly (left-only) -----------------
    Minv_qp = fit_leftOnly_quadPoly(M_qp.Zs, M_qp.ns, sgrid, Minv);

    % ----------------- 8) Fit RM1 and RM2 kernels to quadPoly -----------------
    [Zs_fit, ns_fit] = mergeBasis(FG1_qp.Zs, FG1_qp.ns, FG2_qp.Zs, FG2_qp.ns);
    [Zt_fit, nt_fit] = mergeBasis(FG1_qp.Zt, FG1_qp.nt, FG2_qp.Zt, FG2_qp.nt);
    
    RM1_qp = fit_kernel_from_TS(Zs_fit, ns_fit, Zt_fit, nt_fit, sgrid, T1, S, true);
    RM2_qp = fit_kernel_from_TS(Zs_fit, ns_fit, Zt_fit, nt_fit, sgrid, T2, S, false);

    % ----------------- Diagnostics -----------------
    info.N = N;
    info.h = h;
    info.rankFG1 = n1;
    info.rankFG2 = n2;
    info.n = n;
    info.Lest = Lest;
    info.infoFG1 = infoFG1;
    info.infoFG2 = infoFG2;

    function Ai = A_at(ii)
        if opts.storeABC
            Ai = A(:,:,ii);
        else
            Ai = B(:,:,ii) * C(:,:,ii);
        end
    end
end


% ======================================================================
% Fit kernel quadPoly from T (m×n×N) and S (n×m×N) without storing eta(:,:,i,j)
% ======================================================================
function Q = fit_kernel_from_TS(Zs_out, ns_out, Zt_in, nt_in, grid, T, S, isLower)
% Fit kernel using two-sided least squares:
%   Y ≈ Phi_out * C * Phi_in'
% via QR: Phi_out=Qo*Ro, Phi_in=Qi*Ri, so
%   C = Ro \ (Qo' * Y * Qi) / (Ri')
%
% This avoids kron(Phi_in,Phi_out) and avoids build_Phi().

    N = numel(grid);
    m = size(T,1);
    n = size(T,2);

    % Build basis matrices efficiently (replaces build_Phi)
    Phi_out = basis_matrix_cached(Zs_out, grid);  % N×d_out
    Phi_in  = basis_matrix_cached(Zt_in,  grid);  % N×d_in
    d_out = size(Phi_out,2);
    d_in  = size(Phi_in,2);

    % Thin QR factorizations (reuse for all (p,q))
    [Qo, Ro] = qr(Phi_out, 0);
    [Qi, Ri] = qr(Phi_in,  0);

    % Triangle mask
    [Iidx, Jidx] = ndgrid(1:N, 1:N);
    if isLower
        mask = (Jidx <= Iidx);
    else
        mask = (Jidx > Iidx);
    end

    

    Isp = []; Jsp = []; Vsp = [];

    for p = 1:m
        % A_p(i,:) = T(p,:,i)
        Ap = zeros(N, n);
        for i = 1:N
            Ap(i,:) = T(p,:,i);
        end

        for q = 1:m
            % S_q(:,j) = S(:,q,j)
            Sq = zeros(n, N);
            for j = 1:N
                Sq(:,j) = S(:,q,j);
            end

            % Build Y = Ap * Sq and apply triangle mask
            Y = Ap * Sq;      % N×N
            Y(~mask) = 0;

            % Two-sided LS fit: Cpq = argmin ||Phi_out*C*Phi_in' - Y||_F
            Cpq = Ro \ (Qo' * Y * Qi) / (Ri');

            % Insert Cpq into quadPoly coefficient matrix blocks
            for io = 1:d_out
                row = (p-1)*d_out + io;
                for ii = 1:d_in
                    val = Cpq(io,ii);
                    if val ~= 0
                        col = (q-1)*d_in + ii;
                        Isp(end+1,1) = row; %#ok<AGROW>
                        Jsp(end+1,1) = col; %#ok<AGROW>
                        Vsp(end+1,1) = val; %#ok<AGROW>
                    end
                end
            end
        end
    end
    
    C = sparse(Isp, Jsp, Vsp, m*d_out, m*d_in);
    Q = quadPoly(C,Zs_out,Zt_in,[m,m],ns_out,nt_in);
end

% ======================================================================
% quadPoly evaluation at a scalar point (sval,tval). Use [] for constant side.
% ======================================================================
function F = quadPoly_eval(Q, sval, tval)
    m = Q.dim(1); n = Q.dim(2);

    zs = 1;
    if ~isempty(Q.Zs)
        zs = 1;
        for i = 1:numel(Q.Zs)
            exps = Q.Zs{i}(:);
            if isempty(sval), xi = 1; else, xi = sval; end
            zs = kron(zs, xi.^exps);
        end
    end
    zs = zs(:);

    zt = 1;
    if ~isempty(Q.Zt)
        zt = 1;
        for i = 1:numel(Q.Zt)
            exps = Q.Zt{i}(:);
            if isempty(tval), xi = 1; else, xi = tval; end
            zt = kron(zt, xi.^exps);
        end
    end
    zt = zt(:);

    L = kron(speye(m), zs');   % (m × m*ds)
    R = kron(speye(n), zt);    % (n*dt × n)
    F = full(L * Q.C * R);
end

% ======================================================================
% Fit left-only quadPoly from samples Minv(:,:,i) on grid(i)
% ======================================================================
function Q = fit_leftOnly_quadPoly(Zs, ns, grid, Msamp)
    N = numel(grid);
    m = size(Msamp,1);

    Phi = basis_matrix_cached(Zs, grid);  % N×ds
    ds  = size(Phi,2);

    

    I = []; J = []; V = [];
    for p = 1:m
        for q = 1:m
            y = squeeze(Msamp(p,q,:)); % N×1
            c = Phi \ y;               % ds×1
            for k = 1:ds
                val = c(k);
                if val ~= 0
                    I(end+1,1) = (p-1)*ds + k; %#ok<AGROW>
                    J(end+1,1) = q;           %#ok<AGROW>
                    V(end+1,1) = val;         %#ok<AGROW>
                end
            end
        end
    end
    C = sparse(I,J,V,m*ds,m);
    Q = quadPoly(C,Zs,{0},[m,m],ns,{});
end
function Phi = basis_matrix_cached(Zcell, x)
% Phi(i,:) = Z(x_i)^T, where Z = kron over variables of x.^exps.
% Efficient replacement for build_Phi().

    N = numel(x);
    if isempty(Zcell)
        Phi = ones(N,1);
        return;
    end

    Phi = ones(N,1); % start with scalar basis
    for k = 1:numel(Zcell)
        exps = Zcell{k}(:).';     % 1×nk
        Vk = x(:) .^ exps;        % N×nk
        Phi = rowwise_kron(Phi, Vk);
    end
end

function Z = rowwise_kron(A, B)
% Row-wise Kronecker product:
% If A is N×p and B is N×q, returns Z (N×(p*q)) with
% Z(i,:) = kron(A(i,:), B(i,:)).

    [N,p] = size(A);
    [N2,q] = size(B);
    if N ~= N2
        error('rowwise_kron: row sizes do not match.');
    end

    % Compute all pairwise products per row, then reshape.
    % Result: N×p×q -> N×(p*q)
    Z = reshape(permute(reshape(A, N, p, 1) .* reshape(B, N, 1, q), [1 2 3]), N, p*q);
end



% ======================================================================
% Basis size product
% ======================================================================
function d = quadPoly_basisSize(Zcell)
    if isempty(Zcell), d = 1; return; end
    d = 1;
    for i = 1:numel(Zcell)
        d = d * numel(Zcell{i});
    end
end

% ======================================================================
% Zero kernel quadPoly like another kernel (same bases/dims)
% ======================================================================
function Z = zeroKernelLike(Klike)
    m = Klike.dim(1);
    d_out = quadPoly_basisSize(Klike.Zs);
    d_in  = quadPoly_basisSize(Klike.Zt);

    Z = quadPoly();
    Z.Zs = Klike.Zs; Z.ns = Klike.ns;
    Z.Zt = Klike.Zt; Z.nt = Klike.nt;
    Z.dim = [m m];
    Z.C = sparse(m*d_out, m*d_in);
end
function [Zc, nc] = mergeBasis(Z1, n1, Z2, n2)
%MERGEBASIS Merge quadPoly basis cells (Z, names) variable-by-variable.
% Assumes both sides use the SAME variable ordering and names when present.
% If one side is "constant-only" ({0} or empty), it will be expanded to match the other.
%
% Inputs:
%   Z1, Z2 : 1×k cell arrays of exponent vectors (each a column)
%   n1, n2 : 1×k cell arrays of variable name strings (can be empty)
%
% Outputs:
%   Zc : merged exponent cell array
%   nc : merged variable name cell array

    % Treat empty as constant-only
    if isempty(Z1), Z1 = {0}; end
    if isempty(Z2), Z2 = {0}; end
    if isempty(n1), n1 = repmat({''}, 1, numel(Z1)); end
    if isempty(n2), n2 = repmat({''}, 1, numel(Z2)); end

    % If one is constant-only and the other has multiple variables, expand constant side
    if numel(Z1) == 1 && isConstantBasisCell(Z1) && numel(Z2) > 1
        Z1 = repmat({0}, 1, numel(Z2));
        n1 = repmat({''}, 1, numel(Z2));
    elseif numel(Z2) == 1 && isConstantBasisCell(Z2) && numel(Z1) > 1
        Z2 = repmat({0}, 1, numel(Z1));
        n2 = repmat({''}, 1, numel(Z1));
    end

    if numel(Z1) ~= numel(Z2)
        error('mergeBasis: incompatible number of variables (%d vs %d).', numel(Z1), numel(Z2));
    end

    k = numel(Z1);
    Zc = cell(1,k);
    nc = cell(1,k);

    for i = 1:k
        e1 = Z1{i}(:);
        e2 = Z2{i}(:);

        % Merge and sort unique degrees
        Zc{i} = unique([e1; e2], 'sorted');

        % Pick a name (prefer non-empty)
        name = '';
        if i <= numel(n1) && ~isempty(n1{i}), name = n1{i}; end
        if isempty(name) && i <= numel(n2) && ~isempty(n2{i}), name = n2{i}; end
        nc{i} = name;
    end
end

function tf = isConstantBasisCell(Z)
% True if Z represents constant-only basis (either {0} or a single 0 entry).
    if isempty(Z)
        tf = true;
        return;
    end
    if numel(Z) ~= 1
        tf = false;
        return;
    end
    v = Z{1};
    tf = ~isempty(v) && all(v(:) == 0);
end
