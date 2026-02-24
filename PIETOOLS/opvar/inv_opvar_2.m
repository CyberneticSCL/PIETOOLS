function [Pinv, info] = inv_opvar_2(P, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pinv] = inv_opvar(Pop, opts) computes the inverse operator of the
% operator
% P = [P, Q1; Q2, {R0, R1, R2}];
% INPUT
%   Pop: positive definite opvar to invert
%   opts: inverse options; 
%       - opts.N: number of points in [a,b] to approximate the inverse
%       - opts.mulDeg: maximum degree of multiplier in the inverse
%       - opts.kerDeg: maximum degree of semiseperable kernel in the
%       inverse 
% 
% OUTPUT 
%   Pinv: inverse opvar object. Inverse opvar is a numerical inversion and
%   should be used with care and reservations.
%   info: useful info on conditioning and accuracy of the inverse
% 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - inv_opvar
%
% Copyright (C)2019  M. Peet, S. Shivakumar
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding SS - 2/24/2026
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(opts)
    opts = struct();
end

if ~isa(P,'opvar') || any(P.dim(:,1)~=P.dim(:,2))
    error('Only opvar objects with equal input/output dimensions can be inverted.');
end

% --- options ---
if ~isfield(opts,'N'),        opts.N = 100; end
if ~isfield(opts,'mulDeg'),  opts.mulDeg = 6; end
if ~isfield(opts,'kerDeg'), opts.kerDeg = [6 6]; end

info = struct();

a = P.I(1); b = P.I(2);                % extract the domain and varnames
var1 = P.var1; var2 = P.var2;

% ----- Pure finite-dimensional block (no distributed part) -----
if all(P.dim(2,:)==0)
    if isa(P.P,'polynomial'), P.P = double(P.P); end
    Pinv = P;
    Pinv.P = inv(P.P);
    if nargout>1
        R = struct('cond', cond(P.P));
        info.R = R; 
    end
    return
end

% ----- Pure integral-part operator (no finite-dimensional blocks) -----
if all(P.dim(1,:)==0)

    R0 = P.R.R0;  R1 = P.R.R1;  R2 = P.R.R2;
    
    % Grid for RK4 on [a,b]   
    N = opts.N;
    t = linspace(a,b,N);
    h = t(2)-t(1);
    tmid = t(1:end-1) + 0.5*h;

    % Matrix size
    m = size(R0,1);

    % ---- Evaluate and invert R0 at nodes and midpoints ----
    R0inv_nodes = zeros(m,m,N);
    R0inv_mid   = zeros(m,m,N-1);
    condR0_nodes = zeros(N,1);
    condR0_mid   = zeros(N-1,1);

    for i=1:N
        R0i = double(subs(R0, var1, t(i)));
        condR0_nodes(i) = cond(R0i);
        R0inv_nodes(:,:,i) = inv_via_lu(R0i);
    end
    for i=1:N-1
        R0i = double(subs(R0, var1, tmid(i)));
        condR0_mid(i) = cond(R0i);
        R0inv_mid(:,:,i) = inv_via_lu(R0i);
    end
    
    % ---- Build separable factors for R1 and R2 (polynomial-term stacking) ----
    % Use the polynomial storage to
    % extract monomials once, then evaluate by powers.
    [F1_fun, G1_fun, n1] = factorize_kernel(R1, m, var1.varname, var2.varname);   % splits R1(var1,var2) = F1(var1)*G1(var2), n1 is number of columns in F1
    [F2_fun, G2_fun, n2] = factorize_kernel(R2, m, var1.varname, var2.varname);   % splits R2(var1,var2) = F2(var1)*G2(var2), n2 is number of columns in F2
    n = n1 + n2;

    % ---- Build A(t) = B(t)C(t) for operator I - K, where K uses -R0^{-1}R1/R2 ----
    % For RK4 we evaluate A at nodes and midpoints by index (fast).
    C_node = cell(N,1);
    B_node = cell(N,1);
    A_node = cell(N,1);

    for i=1:N
        C_node{i} = -R0inv_nodes(:,:,i)*[F1_fun(t(i)), F2_fun(t(i))];  % m x n,  C = -R0^{-1}*[F1, F2]
        B_node{i} = [G1_fun(t(i)); -G2_fun(t(i))]; % n x m, B = [G1; -G2]
        A_node{i} = B_node{i}*C_node{i};       % n x n
    end
    
    % repeat for mid points of the stencil
    C_mid = cell(N-1,1);
    B_mid = cell(N-1,1);
    A_mid = cell(N-1,1);

    for i=1:N-1
        C_mid{i} = -R0inv_mid(:,:,i)*[F1_fun(tmid(i)), F2_fun(tmid(i))];
        B_mid{i} = [G1_fun(tmid(i)); -G2_fun(tmid(i))];
        A_mid{i} = B_mid{i} * C_mid{i};
    end

    % ---- RK4 propagate U and V=U^{-1} using precomputed A at nodes/mids ----
    % we propagate by solving d/dt U(t) = A(t)U(t); U(0) = I
    %                         d/dt V(t) = -V(t)A(t); V(0) = I
    [U, V] = rk4_UV(A_node, A_mid, t);

    % ---- Indicator block & P matrix ----
    % split U(b) = [U11(b), U12(b); U21(b), U22(b)]
    Ub  = U(:,:,end);
    U21 = Ub(n1+1:end,1:n1);
    U22 = Ub(n1+1:end,n1+1:end);

    
    if rcond(U22) < 1e-12
        warning('U22(b) ill-conditioned; inverse may be unstable.');
    end

    X = U22 \ U21;  % n2 x n1
    Pmat = [zeros(n1,n1), zeros(n1,n2);
            X,            eye(n2)];
    ImP  = eye(n) - Pmat;

    % ---- Precompute M_i = C(t_i)U(t_i), N_j = V(t_j)B(t_j) ----
    % only at nodes, since midpoints were used purely for rk4 integration
    M = zeros(m,n,N);
    Nmat = zeros(n,m,N);

    for i=1:N
        M(:,:,i)    = C_node{i}*U(:,:,i);
        Nmat(:,:,i) = V(:,:,i)*B_node{i};
    end

    % ---- Build S0,S1,S2 numerically at nodes ----
    S0_num = R0inv_nodes;  % m x m x N, S0 = R0^{-1}

    % S1_num, S2_num: m x m x N x N
    S1_num = zeros(m,m,N,N);   % S1 = C*U*(I-P)*V*B*R0^{-1}
    S2_num = zeros(m,m,N,N);   % S2 = -C*U*(P)*V*B*R0^{-1}
    
    if exist('pagemtimes','file') == 2
        % Vectorized over j for each i
        for i=1:N
            Mi = M(:,:,i);

            % lower: gamma = Mi*ImP*Nmat(:,:,j), j<=i
            Glo = pagemtimes(Mi*ImP, Nmat);            % m x m x N
            % upper: gamma = -Mi*Pmat*Nmat(:,:,j), j>i
            Gup = pagemtimes(-Mi*Pmat, Nmat);          % m x m x N

            % Multiply by R0inv(s)^{-1} on each page j:
            Slo = pagemtimes(Glo, R0inv_nodes);        % m x m x N
            Sup = pagemtimes(Gup, R0inv_nodes);        % m x m x N

            % Assign triangular parts
            if i >= 1
                S1_num(:,:,i,1:i) = reshape(Slo(:,:,1:i), m,m,1,i);
            end
            if i < N
                S2_num(:,:,i,i+1:N) = reshape(Sup(:,:,i+1:N), m,m,1,N-i);
            end
        end
    else
        % Fallback loops (slower)
        for i=1:N
            for j=1:N
                if t(j) <= t(i)
                    gamma_ij = M(:,:,i)*ImP*Nmat(:,:,j);
                    S1_num(:,:,i,j) = gamma_ij*R0inv_nodes(:,:,j);
                else
                    gamma_ij = -M(:,:,i)*Pmat*Nmat(:,:,j);
                    S2_num(:,:,i,j) = gamma_ij * R0inv_nodes(:,:,j);
                end
            end
        end
    end

    % ---- Fit polynomials at the end ----
    [S0_poly, fitS0] = fitpoly_1D_cheb(t, S0_num, opts.mulDeg, var1.varname, [a b]); % convert (m,m,N) matrix of points to (m,m) polynomial in var1
    [S1_poly, fitS1] = fitpoly_2D_cheb(t, t, S1_num, opts.kerDeg, var1.varname, var2.varname, 'lower', [a b], [a b]); % convert (m,m,N,N) matrix of points to (m,m) polynomial in (var1,var2)
    [S2_poly, fitS2] = fitpoly_2D_cheb(t, t, S2_num, opts.kerDeg, var1.varname, var2.varname, 'upper', [a b], [a b]); % convert (m,m,N,N) matrix of points to (m,m) polynomial in (var1,var2)
    
    Pinv = P;
    Pinv.R.R0 = S0_poly;
    Pinv.R.R1 = S1_poly;
    Pinv.R.R2 = S2_poly;
    
    if nargout>1
        infoL2.condR0.max = max([condR0_nodes(:);condR0_mid(:)]);
        infoL2.condR0.min = min([condR0_nodes(:);condR0_mid(:)]);
        infoL2.fitS0 = fitS0;
        infoL2.fitS1 = fitS1;
        infoL2.fitS2 = fitS2;
        infoL2.condU22 = cond(U22);
        info.L2 = infoL2;
    end

    return
end

% ----- Full 2x2 block opvar case (your Schur complement logic) -----
opvar A B C D;
A.I = P.I; A.var1 = P.var1; A.var2 = P.var2;
B.I = P.I; B.var1 = P.var1; B.var2 = P.var2;
C.I = P.I; C.var1 = P.var1; C.var2 = P.var2;
D.I = P.I; D.var1 = P.var1; D.var2 = P.var2;

A.P  = P.P;
B.Q1 = P.Q1;
C.Q2 = P.Q2;
D.R  = P.R;

[Ainv, infoFinite] = inv_opvar_2(A, opts);
[TB, infoInfinite]   = inv_opvar_2(D - C*Ainv*B, opts);

Pinv = [Ainv + Ainv*B*TB*C*Ainv,  -Ainv*B*TB;
        -TB*C*Ainv,                TB];

if nargout>1
    info = [infoFinite, infoInfinite];
end
end

function X = inv_via_lu(A)
% Return A^{-1} via LU solves (no explicit inv(A) call).
[L,U,P] = lu(A);
X = U \ (L \ P);
end

function [U,V] = rk4_UV(A_node, A_mid, tgrid)
% A_node{k} = A(t_k), k=1..N
% A_mid{k}  = A(t_k + h/2), k=1..N-1
n = size(A_node{1},1);
N = numel(tgrid);

U = zeros(n,n,N);
V = zeros(n,n,N);
U(:,:,1) = eye(n);
V(:,:,1) = eye(n);

h = tgrid(2) - tgrid(1);

for k=1:N-1
    Uk = U(:,:,k);
    Vk = V(:,:,k);

    A1 = A_node{k};
    A2 = A_mid{k};
    A4 = A_node{k+1};

    % U' = A U
    k1U = A1*Uk;
    k2U = A2*(Uk + 0.5*h*k1U);
    k3U = A2*(Uk + 0.5*h*k2U);
    k4U = A4*(Uk + h*k3U);

    % V' = -V A
    k1V = -Vk*A1;
    k2V = -(Vk + 0.5*h*k1V)*A2;
    k3V = -(Vk + 0.5*h*k2V)*A2;
    k4V = -(Vk + h*k3V)*A4;

    U(:,:,k+1) = Uk + (h/6)*(k1U + 2*k2U + 2*k3U + k4U);
    V(:,:,k+1) = Vk + (h/6)*(k1V + 2*k2V + 2*k3V + k4V);
end
end

function [P, info] = fitpoly_1D_cheb(t, Msamp, degT, varname, interval)
%FITPOLY_1D_CHEB  Fit matrix-valued polynomial using Chebyshev basis (stable),
% then convert to monomial-basis 'polynomial' object.
%
%   [P, info] = fitpoly_1D_cheb(t, Msamp, degT, varname, interval)
%
% Inputs
%   t        : 1xN or Nx1 sample points (assumed in [a,b])
%   Msamp    : m x n x N samples
%   degT     : max degree (>=0)
%   varname  : variable name (string)
%   interval : [a b] (optional). If omitted uses [min(t) max(t)].
%
% Outputs
%   P        : polynomial object in monomial basis (as required by your class)
%   info     : struct with fields:
%                .rms_residual, .rel_rms_residual, .cond

N = numel(t);
[m,n,NN] = size(Msamp);
if NN ~= N, error('Msamp third dimension must match numel(t).'); end

degT = round(degT);
if degT < 0, error('degT must be nonnegative.'); end

if nargin < 5 || isempty(interval)
    a = min(t); b = max(t);
else
    a = interval(1); b = interval(2);
end
if ~(b > a), error('interval must satisfy b>a.'); end

% Map t in [a,b] to x in [-1,1]
x = (2*(t - a)/(b - a)) - 1;

% Chebyshev design matrix: T_k(x), k=0..degT
V = chebV_1D(x, degT);                 % N x (degT+1)

% Stack samples as N x (m*n)
Y = reshape(permute(Msamp, [3 1 2]), N, m*n);

% Solve least squares in Cheb basis
Ccheb = V \ Y;                         % (degT+1) x (m*n)

% Residual stats (effective)
R = V*Ccheb - Y;
rms = sqrt(mean(R(:).^2));
relrms = rms / max(1e-14, sqrt(mean(Y(:).^2)));

info = struct();
info.rms_residual = rms;
info.rel_rms_residual = relrms;
info.cond = cond(V);

% Convert Cheb coefficients to monomial coefficients in t for each matrix entry
Cmono = chebcoef_to_monomial_1D(Ccheb, degT, a, b);  % (degT+1) x (m*n)

degmat = sparse((0:degT).');
coef   = sparse(Cmono);
P = polynomial(coef, degmat, varname, [m n]);
end


function [P, info] = fitpoly_2D_cheb(tgrid, sgrid, Msamp, degTS, tname, sname, region, interval_t, interval_s)
%FITPOLY_2D_CHEB  Fit matrix-valued polynomial using tensor Chebyshev basis on a region,
% then convert to monomial-basis 'polynomial' object.
%
%   [P, info] = fitpoly_2D_cheb(tgrid, sgrid, Msamp, degTS, tname, sname, region, interval_t, interval_s)
%
% Inputs
%   tgrid, sgrid : sample points (vectors)
%   Msamp        : m x n x Nt x Ns samples
%   degTS        : [degT degS]
%   tname,sname  : variable names (strings)
%   region       : 'all'|'lower'|'upper'
%   interval_t   : [a b] for t mapping to [-1,1] (optional)
%   interval_s   : [c d] for s mapping to [-1,1] (optional)
%
% Outputs
%   P            : polynomial object in monomial basis in variables {tname;sname}
%   info         : residual/conditioning stats

Nt = numel(tgrid); Ns = numel(sgrid);

[m,n,NNt,NNs] = size(Msamp);
if NNt~=Nt || NNs~=Ns
    error('Msamp must be m x n x Nt x Ns matching tgrid and sgrid.');
end

degT = round(degTS(1));
degS = round(degTS(2));
if degT < 0 || degS < 0, error('degTS must be nonnegative.'); end

if nargin < 7 || isempty(region), region = 'all'; end
region = lower(region);

if nargin < 8 || isempty(interval_t)
    at = min(tgrid); bt = max(tgrid);
else
    at = interval_t(1); bt = interval_t(2);
end
if nargin < 9 || isempty(interval_s)
    as = min(sgrid); bs = max(sgrid);
else
    as = interval_s(1); bs = interval_s(2);
end
if ~(bt > at) || ~(bs > as), error('intervals must satisfy upper>lower.'); end

% Build all pairs and apply region mask
[Tmat,Smat] = ndgrid(tgrid, sgrid);   % Nt x Ns
switch region
    case 'all'
        mask = true(Nt,Ns);
    case 'lower'
        mask = (Smat <= Tmat);
    case 'upper'
        mask = (Smat > Tmat);
    otherwise
        error('region must be all|lower|upper');
end

tt = Tmat(mask);
ss = Smat(mask);
Mpts = numel(tt);

% Map to [-1,1]
xt = (2*(tt - at)/(bt - at)) - 1;
xs = (2*(ss - as)/(bs - as)) - 1;

% 1D Cheb matrices at selected points
Vt = chebV_1D(xt, degT);              % Mpts x (degT+1)
Vs = chebV_1D(xs, degS);              % Mpts x (degS+1)

% Build tensor-product design matrix Phi where column corresponds to (p,q):
% Vt: Mpts x (degT+1), Vs: Mpts x (degS+1)
% Phi: Mpts x ((degT+1)*(degS+1))
% Phi(:, idx(p,q)) = T_p(xt) .* T_q(xs)
Phi = kron(Vt, ones(1,degS+1)).*kron(ones(1,degT+1), Vs);

% Stack samples into Y (Mpts x (m*n))
Yfull = reshape(permute(Msamp, [3 4 1 2]), Nt*Ns, m*n);
Y = Yfull(mask(:), :);

% Solve least squares in Cheb-tensor basis
Ccheb = Phi\Y;                       % Tterms x (m*n)

% Residual stats
R = Phi*Ccheb - Y;
rms = sqrt(mean(R(:).^2));
relrms = rms/max(1e-14, sqrt(mean(Y(:).^2)));

info = struct();
info.rms_residual = rms;
info.rel_rms_residual = relrms;
info.cond = cond(Phi);


% Convert Cheb-tensor coefficients to monomial in (t,s)
Cmono = chebcoef_to_monomial_2D(Ccheb, degT, degS, at, bt, as, bs);  % Tterms x (m*n)

% Build polynomial (monomial basis degrees)
[pdeg,qdeg] = ndgrid(0:degT, 0:degS);
pdeg = pdeg(:); qdeg = qdeg(:);
degmat = sparse([pdeg,qdeg]);
coef = sparse(Cmono);

P = polynomial(coef, degmat, [tname; sname], [m n]);
end


% ------------------------- helpers -------------------------

function V = chebV_1D(x, deg)
% Return matrix V where V(:,k+1) = T_k(x), k=0..deg
x = x(:);
N = numel(x);
V = zeros(N, deg+1);
V(:,1) = 1;
if deg >= 1
    V(:,2) = x;
end
for k = 2:deg
    V(:,k+1) = 2*x.*V(:,k) - V(:,k-1);
end
end

function Cmono = chebcoef_to_monomial_1D(Ccheb, deg, a, b)
% Convert Cheb coeffs in x to monomial coeffs in t.
% x = alpha*t + beta, alpha = 2/(b-a), beta = -(a+b)/(b-a).
%
% Input:
%   Ccheb : (deg+1) x (m*n), row k+1 corresponds to T_k(x)
% Output:
%   Cmono : (deg+1) x (m*n), row r+1 corresponds to t^r
alpha = 2/(b-a);
beta  = -(a+b)/(b-a);

% Build conversion matrix M such that:
%   [t^0; t^1; ... t^deg] coefficients = M * [T_0; ... T_deg] coefficients
% i.e., for scalar: sum c_k T_k(alpha t + beta) = sum d_r t^r, with d = M*c
M = cheb_to_monomial_matrix_1D(deg, alpha, beta);   % (deg+1) x (deg+1)
Cmono = M * Ccheb;                                  % (deg+1) x (m*n)
end


function Cmono = chebcoef_to_monomial_2D(Ccheb, degT, degS, at, bt, as, bs)
% Convert tensor Cheb coeffs to monomial coeffs in (t,s).
% Ordering of Ccheb rows matches:
%   col = 1;
%   for p=0:degT
%     for q=0:degS
%        row corresponds to T_p(xt)*T_q(xs)
% xt = alpha_t*t + beta_t, xs = alpha_s*s + beta_s

alpha_t = 2/(bt-at); beta_t = -(at+bt)/(bt-at);
alpha_s = 2/(bs-as); beta_s = -(as+bs)/(bs-as);

Mt = cheb_to_monomial_matrix_1D(degT, alpha_t, beta_t); % (degT+1)x(degT+1)
Ms = cheb_to_monomial_matrix_1D(degS, alpha_s, beta_s); % (degS+1)x(degS+1)

% We want to map coefficients c_{p,q} in Cheb basis to monomial d_{i,j}.
% In matrix form for scalar:
%   D = Mt * C * Ms^T
% where C is (degT+1)x(degS+1) Cheb coefficients, D is monomial coeffs.
%
% Here we have (m*n) columns, each column is a scalar surface coefficient set.
Tterms = (degT+1)*(degS+1);
mn = size(Ccheb,2);
Cmono = zeros(Tterms, mn);

for col = 1:mn
    Cpq = reshape(Ccheb(:,col), [degS+1, degT+1]).';   % (degT+1)x(degS+1)
    Dij = Mt * Cpq * (Ms.');                           % (degT+1)x(degS+1)
    Cmono(:,col) = Dij(:);                             % match ndgrid p,q ordering after vectorize
end
end


function M = cheb_to_monomial_matrix_1D(deg, alpha, beta)
% Build M (deg+1 x deg+1) such that:
%   sum_{k=0}^deg c_k T_k(alpha t + beta) = sum_{r=0}^deg d_r t^r,
% with d = M*c.
%
% Coeff vectors are ascending powers of t: [c0 c1 ...].
Tpoly = cell(deg+1,1);

% T0 = 1
Tpoly{1} = 1;                 % row vector

if deg >= 1
    % T1(x)=x, x = beta + alpha t
    Tpoly{2} = [beta, alpha]; % row vector
end

for k = 2:deg
    % T_k(x) = 2 x T_{k-1}(x) - T_{k-2}(x)
    p = Tpoly{k};     % coeffs for T_{k-1}
    q = Tpoly{k-1};   % coeffs for T_{k-2}

    p = p(:).';       % FORCE row vector
    q = q(:).';       % FORCE row vector

    % x*p where x = beta + alpha t
    xp = conv(p, [beta, alpha]);      % row vector

    qpad = pad_to_row(q, numel(xp));  % row vector same length as xp
    Tk = 2*xp - qpad;                 % row vector
    Tpoly{k+1} = Tk;
end

% Assemble M: column (k+1) = coeffs of T_k(alpha t + beta) padded/truncated to deg
M = zeros(deg+1, deg+1);
for k = 0:deg
    coeffs = Tpoly{k+1};
    coeffs = coeffs(:).';                 % row
    coeffs = pad_to_row(coeffs, deg+1);   % row length deg+1
    M(:,k+1) = coeffs(:);                 % column
end
end

function v = pad_to_row(v, L)
% v row vector, pad or truncate to length L
v = v(:).';
nv = numel(v);
if nv < L
    v = [v, zeros(1, L-nv)];
elseif nv > L
    v = v(1:L);
end
end

function [F_fun, G_fun, nFG] = factorize_kernel(P, m, tname, sname)
terms = poly_terms_2D(P, tname, sname);
T = numel(terms);
nFG = m*T;

F_fun = @(tt) F_eval(tt, terms, m);
G_fun = @(ss) G_eval(ss, terms, m);
end

%cheaper subs function for polynomial class objects
function F = F_eval(tt, terms, m)
T = numel(terms);
F = zeros(m, m*T);
for k=1:T
    cols = (k-1)*m + (1:m);
    F(:,cols) = (tt.^terms(k).p) * terms(k).C;
end
end

function G = G_eval(ss, terms, m)
T = numel(terms);
G = zeros(m*T, m);
Im = eye(m);
for k=1:T
    rows = (k-1)*m + (1:m);
    G(rows,:) = (ss.^terms(k).q) * Im;
end
end

function terms = poly_terms_2D(P, tname, sname)
V = numel(P.varname);
if V ~= 2, error('Expected a 2-variable polynomial.'); end

idxT = find(strcmp(P.varname, tname), 1);
idxS = find(strcmp(P.varname, sname), 1);
if isempty(idxT) || isempty(idxS)
    error('Variables not found in P.varname.');
end

T = size(P.degmat,1);
m = P.matdim(1); n = P.matdim(2);
if m ~= n
    error('This separable builder assumes square m x m matrices.');
end

D = full(P.degmat);
Cmat = full(P.coefficient);  % T x (m*m)

terms = repmat(struct('p',0,'q',0,'C',zeros(m,m)), T, 1);
for k=1:T
    terms(k).p = D(k,idxT);
    terms(k).q = D(k,idxS);
    terms(k).C = reshape(Cmat(k,:), [m,m]); % column-major matches class comment
end
end

