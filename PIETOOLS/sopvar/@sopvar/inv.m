function [Rinv,info] = inv(R,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [RINV,INFO] = INV(R[,OPTS]) the inverse of a square 1-D 'sopvar' block
%   (R x)(s) = R0(s) x(s) + int_a^s R1(s,t) x(t) dt + int_s^b R2(s,t) x(t) dt
% on L2^m[a,b], as a 'sopvar' with polynomial parameters of one degree d in
% s and in t. The block stores (sopvar document Sec. 4, 1-D)
%   R_i(s,t) = (I_m kron ZL(s)') C_i (I_m kron ZR(t)),   C_i = R.params{i+1},
% so each kernel is already in the separated form Lemma 16 below needs:
%   C_i = U_i V_i  (SVD, numerical rank r_i)  ==>  R_i = F_i(s) G_i(t),
%   F_i(s) = (I kron ZL(s)') U_i   (m x r_i),   G_i(t) = V_i (I kron ZR(t)).
%
% METHOD (Shivakumar, Das, Peet, arXiv 2208.13104, Sec. VII, Lemma 16 and
% Cor. 17, Gohberg-Krein). With H_i = R0^{-1} R_i = -F_i^h G_i, F_i^h =
% -R0^{-1} F_i, set C(s) = [F_1^h F_2^h] (m x n), B(s) = [G_1; -G_2] (n x m),
% n = r_1 + r_2, and solve on [a,b]
%   U' = B C U, U(a) = I;      V' = -V B C, V(a) = I   (V = U^{-1}).
% Partition U(b) with U22 of size r_2; R is invertible iff U22 is, and
%   Pm = [0 0; U22^{-1} U21 I],
%   R^{-1} = {R0^{-1},  C(s) U(s) (I - Pm) V(t) B(t) R0(t)^{-1},
%                      -C(s) U(s) Pm V(t) B(t) R0(t)^{-1}}.
% U, V are propagated by RK4 on OPTS.N uniform nodes with midpoints. The
% three parameters are evaluated on the node grid (the kernels on the WHOLE
% square from the formulas, not only on their triangle, because the formula
% is the smooth continuation and a polynomial fits it; fitting a kernel that
% is cut at the diagonal does not converge) and fitted by least squares in
% the Chebyshev basis on [a,b], degree d raised from OPTS.deg0 in steps of 2
% until the relative RMS residual of all three fits is below OPTS.tol or d =
% OPTS.degmax; the Chebyshev coefficients are converted to the monomial
% basis ZL = ZR = (0:d)', and the multiplier is placed in the ZR-degree-0
% columns, the canonical multiplier form of the class. R0^{-1} is rational,
% so the inverse is only approximated by a polynomial block; INFO reports
% the fit residuals. This is the route of the stock 'inv_opvar_2' (SS,
% 02/2026), written against the block storage: no polynomial objects, no
% monomial bookkeeping, the kernel split is one SVD per kernel.
%
% INPUTS
% - R:      'sopvar', L2^m[s] -> L2^m[s] with one spatial variable, square;
% - opts:   (optional) struct, fields
%     N        nodes of the RK4 grid, default 101 (the midpoints are added);
%     tol      target relative RMS fit residual, default 1e-8;
%     deg0     first fit degree, default 4;
%     degmax   last fit degree tried, default 16;
%     deg      a fixed fit degree (overrides deg0, degmax, tol);
%     tol_rank relative SVD cut for the kernel rank, default 1e-13.
% OUTPUTS
% - Rinv:   'sopvar', the inverse block on ZL = ZR = (0:d)';
% - info:   struct: N, h, r1, r2 (kernel ranks), n, rcondU22, condR0 (max over
%           the grid), d, relrms (1 x 3, the fit residuals of R0^{-1}, R1^h,
%           R2^h), and note.
%
% NOTES
% Cost: O(N n^3) for the RK4 sweep, O(N^2 m^2 n) for the kernel values,
% O(N^2 (d+1)^2 m^2) per fit degree. No decision variables are involved; the
% block is fixed, so nothing here scales with q.
% A singular R0(s) at a grid point, or a near-singular U22(b), is reported
% through INFO and a warning; the formulas then give a large residual rather
% than an error.
%
% See also COPVAR/INV, INV_OPVAR_2, INV_OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - inv (sopvar)
%
% Copyright (C) 2026 PIETOOLS Team
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
% Initial coding MMP, 10/09/2026. Container counterpart of @opvar/inv for
%                the synthesis reconstruction K = Z P^{-1}, L = P^{-1} Z on
%                the container path (getController_sop, getObserver_sop).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    opts = struct();    end
o = struct('N',101,'tol',1e-8,'deg0',4,'degmax',16,'deg',[],'tol_rank',1e-13);
fn = fieldnames(opts);
for k = 1:numel(fn),    o.(fn{k}) = opts.(fn{k});   end
if ~isempty(o.deg),     o.deg0 = o.deg;     o.degmax = o.deg;   end

% ---- the block must be 1-D, L2 -> L2 in one variable, and square
if numel(R.vars.in)~=1 || numel(R.vars.out)~=1 || ~strcmp(R.vars.in{1},R.vars.out{1})
    error('sopvar:inv:dim','INV is implemented for 1-D blocks L2^m[s] -> L2^m[s]; use COPVAR/INV for a grid with an R^n space.');
end
m = R.dims(1);
if R.dims(2)~=m,    error('sopvar:inv:square','The block must be square (dims %d x %d).',R.dims(1),R.dims(2));    end
a = R.dom.in(1,1);  b = R.dom.in(1,2);
zL = R.ZL{1}(:);    zR = R.ZR{1}(:);     nL = numel(zL);  nR = numel(zR);
C0 = zero_or(R.params{1},m*nL,m*nR);
C1 = zero_or(R.params{2},m*nL,m*nR);
C2 = zero_or(R.params{3},m*nL,m*nR);
Im = speye(m);

% ---- kernel factors from the coefficient matrices (one SVD each)
[U1,V1] = factor_C(C1,o.tol_rank);  r1 = size(U1,2);
[U2,V2] = factor_C(C2,o.tol_rank);  r2 = size(U2,2);
n = r1 + r2;

% ---- grid, R0^{-1} and A = B C on nodes and midpoints
N = o.N;    sg = linspace(a,b,N);   h = sg(2) - sg(1);    smid = sg(1:end-1) + h/2;
R0inv = zeros(m,m,N);   condR0 = zeros(N,1);
Cn = cell(N,1);     Bn = cell(N,1);     An = cell(N,1);     Am = cell(N-1,1);
for i = 1:N
    [Li,Ri] = bases(sg(i));
    R0i = full(Li*C0*Ri);   condR0(i) = cond(R0i);  R0inv(:,:,i) = R0i\eye(m);
    [Cn{i},Bn{i}] = CB(Li,Ri,R0inv(:,:,i));
    An{i} = Bn{i}*Cn{i};
end
for i = 1:N-1
    [Li,Ri] = bases(smid(i));
    R0i = full(Li*C0*Ri);
    [Ci,Bi] = CB(Li,Ri,R0i\eye(m));
    Am{i} = Bi*Ci;
end

% ---- RK4: U' = A U, V' = -V A (V = U^{-1}), from U(a) = V(a) = I
U = zeros(n,n,N);   V = zeros(n,n,N);
U(:,:,1) = eye(n);  V(:,:,1) = eye(n);
for k = 1:N-1
    A1 = An{k};     A2 = Am{k};     A4 = An{k+1};
    Uk = U(:,:,k);  Vk = V(:,:,k);
    k1 = A1*Uk;     k2 = A2*(Uk + h/2*k1);  k3 = A2*(Uk + h/2*k2);  k4 = A4*(Uk + h*k3);
    U(:,:,k+1) = Uk + h/6*(k1 + 2*k2 + 2*k3 + k4);
    k1 = -Vk*A1;    k2 = -(Vk + h/2*k1)*A2; k3 = -(Vk + h/2*k2)*A2; k4 = -(Vk + h*k3)*A4;
    V(:,:,k+1) = Vk + h/6*(k1 + 2*k2 + 2*k3 + k4);
end

% ---- Pm from U(b): the upper kernel's block decides invertibility
Ub = U(:,:,N);
U21 = Ub(r1+1:n,1:r1);  U22 = Ub(r1+1:n,r1+1:n);
if r2>0,    rc = rcond(U22);    else,   rc = 1;     end
if rc<1e-14
    warning('sopvar:inv:singular','U22(b) is near singular (rcond %.2g): the block is likely not invertible.',rc);
end
Pm = zeros(n);
if r2>0,    Pm(r1+1:n,1:r1) = U22\U21;  Pm(r1+1:n,r1+1:n) = eye(r2);   end
ImP = eye(n) - Pm;

% ---- parameter values on the node grid: S0_i = R0^{-1}(s_i),
%      S1_ij = M_i (I-Pm) N_j R0inv_j, S2_ij = -M_i Pm N_j R0inv_j, with
%      M_i = C(s_i) U(s_i), N_j = V(s_j) B(s_j); whole square, see header
Mn = zeros(m,n,N);  Nn = zeros(n,m,N);
for i = 1:N
    Mn(:,:,i) = Cn{i}*U(:,:,i);
    Nn(:,:,i) = V(:,:,i)*Bn{i}*R0inv(:,:,i);
end
S1 = zeros(m,m,N,N);    S2 = zeros(m,m,N,N);
if n>0
    for i = 1:N
        G1 = pagemtimes(Mn(:,:,i)*ImP,Nn);          % m x m x N over j
        G2 = pagemtimes(-Mn(:,:,i)*Pm,Nn);
        S1(:,:,i,:) = reshape(G1,m,m,1,N);
        S2(:,:,i,:) = reshape(G2,m,m,1,N);
    end
end

% ---- fit: Chebyshev least squares at one degree d for all three, raised
%      until every relative RMS residual is below tol
x = (2*(sg - a)/(b - a)) - 1;           % nodes mapped to [-1,1]
Y0 = reshape(permute(R0inv,[3 1 2]),N,m*m);             % N x m^2, column-major entries
Y1 = reshape(permute(S1,[3 4 1 2]),N*N,m*m);            % (i outer? no: i inner, j outer) see below
Y2 = reshape(permute(S2,[3 4 1 2]),N*N,m*m);
% permute [3 4 1 2] orders the rows as (i,j) with i fastest; the design
% matrix below is built with the same order (kron over j outer, i inner).
d = o.deg0;
while true
    T1 = chebV(x,d);                                    % N x (d+1), T_p(x_i)
    % the 2-D samples run over (i,j) with i fastest, so Phi(i + (j-1)N,(p,q)) = T_p(x_i) T_q(x_j)
    Phi = design2(T1,d);                                % N^2 x (d+1)^2, row i + (j-1) N, column (p-1)(d+1) + q
    [c0,e0] = lsfit(T1,Y0);
    if n>0
        [c1,e1] = lsfit(Phi,Y1);    [c2,e2] = lsfit(Phi,Y2);
    else
        c1 = zeros((d+1)^2,m*m);    c2 = c1;    e1 = 0;    e2 = 0;
    end
    if max([e0 e1 e2])<=o.tol || d>=o.degmax,   break,  end
    d = min(d+2,o.degmax);
end
% Chebyshev (in x) -> monomial (in s) coefficients: d_r = sum_p M(r,p) c_p
Mc = cheb2mono(d,2/(b-a),-(a+b)/(b-a));                 % (d+1) x (d+1)
C0mono = Mc*c0;                                         % (d+1) x m^2, row r+1 <-> s^r
C1mono = zeros((d+1)^2,m*m);    C2mono = C1mono;
for col = 1:m*m
    Cpq = reshape(c1(:,col),d+1,d+1).';                 % (p,q) with p the s index: rows were (p-1)(d+1)+q
    D = Mc*Cpq*Mc.';                                    % (r,u): s^r t^u
    C1mono(:,col) = reshape(D.',[],1);                  % back to (r-1)(d+1)+u
    Cpq = reshape(c2(:,col),d+1,d+1).';
    D = Mc*Cpq*Mc.';
    C2mono(:,col) = reshape(D.',[],1);
end
% ---- place into the block storage on ZL = ZR = (0:d)'
nZ = d + 1;
C0h = zeros(m*nZ,m*nZ);     C1h = C0h;     C2h = C0h;
for i = 1:m
    for j = 1:m
        col = i + (j-1)*m;                              % column-major entry (i,j)
        C0h((i-1)*nZ+(1:nZ),(j-1)*nZ+1) = C0mono(:,col);        % ZR degree 0 column: canonical multiplier form
        C1h((i-1)*nZ+(1:nZ),(j-1)*nZ+(1:nZ)) = reshape(C1mono(:,col),nZ,nZ).';   % (r,u) -> row r, col u
        C2h((i-1)*nZ+(1:nZ),(j-1)*nZ+(1:nZ)) = reshape(C2mono(:,col),nZ,nZ).';
    end
end
Zh = (0:d)';
Rinv = sopvar({sparse(C0h);sparse(C1h);sparse(C2h)},R.vars,{Zh},{Zh},R.dom,[m m]);

info = struct('N',N,'h',h,'r1',r1,'r2',r2,'n',n,'rcondU22',rc,'condR0',max(condR0),...
              'd',d,'relrms',[e0 e1 e2],'note','Lemma 16 / Cor. 17 of arXiv 2208.13104, RK4 + Chebyshev LS fit');
if max(info.relrms)>o.tol
    info.note = sprintf('%s; fit tolerance %.1e NOT met at degmax %d (relrms %.2e)',info.note,o.tol,o.degmax,max(info.relrms));
end

    function [L,Rr] = bases(s)
        L = kron(Im,sparse(s.^zL)');    Rr = kron(Im,sparse(s.^zR));
    end
    function [Ci,Bi] = CB(L,Rr,R0i)
        % C(s) = -R0(s)^{-1} [F1 F2], B(s) = [G1; -G2]
        F = [full(L*U1), full(L*U2)];
        G = [full(V1*Rr); -full(V2*Rr)];
        Ci = -R0i*F;    Bi = G;
    end
end


function C = zero_or(C,nr,nc)
% [] and scalar 0 are the zero-block shorthand of the class
if isempty(C) || (isscalar(C) && C==0),    C = sparse(nr,nc);   end
if ~isequal(size(C),[nr nc])
    error('sopvar:inv:params','A parameter is %d x %d; the bases give %d x %d.',size(C,1),size(C,2),nr,nc);
end
end

function [U,V] = factor_C(C,tol)
% C = U V with the numerical rank: economy SVD, singular values below tol*max dropped
if nnz(C)==0,   U = zeros(size(C,1),0);  V = zeros(0,size(C,2));  return,  end
[Us,S,Vs] = svd(full(C),'econ');
sv = diag(S);   r = find(sv>=tol*sv(1),1,'last');
U = Us(:,1:r)*diag(sv(1:r));    V = Vs(:,1:r)';
end

function T = chebV(x,d)
% T(:,p+1) = T_p(x), p = 0..d
x = x(:);   T = zeros(numel(x),d+1);    T(:,1) = 1;
if d>=1,    T(:,2) = x;     end
for p = 2:d,    T(:,p+1) = 2*x.*T(:,p) - T(:,p-1);   end
end

function Phi = design2(T1,d)
% rows (i,j) with i fastest, columns (p-1)(d+1)+q: Phi = T_p(x_i) T_q(x_j)
N = size(T1,1);
Phi = zeros(N*N,(d+1)^2);
for q = 1:d+1
    for p = 1:d+1
        Phi(:,(p-1)*(d+1)+q) = reshape(T1(:,p)*T1(:,q)',[],1);   % (i,j) -> i + (j-1) N
    end
end
end

function [c,relrms] = lsfit(Phi,Y)
c = Phi\Y;
res = Phi*c - Y;
relrms = sqrt(mean(res(:).^2))/max(1e-300,sqrt(mean(Y(:).^2)));
if ~any(Y(:)),  relrms = 0;     end
end

function M = cheb2mono(d,alpha,beta)
% sum_p c_p T_p(alpha s + beta) = sum_r (M c)_r s^r ; column p+1 of M holds T_p(alpha s + beta)
M = zeros(d+1,d+1);
Tprev = 1;  M(1,1) = 1;
if d>=1
    Tcur = [beta, alpha];   M(1:2,2) = Tcur(:);
    for p = 2:d
        Tnext = 2*conv(Tcur,[beta, alpha]);
        Tnext(1:numel(Tprev)) = Tnext(1:numel(Tprev)) - Tprev;
        M(1:numel(Tnext),p+1) = Tnext(:);
        Tprev = Tcur;   Tcur = Tnext;
    end
end
end
