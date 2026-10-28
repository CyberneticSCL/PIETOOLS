function [K,info] = getController_direct_sop(P,Z,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [K,INFO] = GETCONTROLLER_DIRECT_SOP(P,Z[,OPTS]) returns the
% state-feedback gain
%   u = K x = Z P^{-1} x,   x = [x1; x2(.)] in R^k x L2^m[a,b],
% without forming P^{-1} as an operator: every matrix the gain needs is an
% integral of grid values of the Gohberg-Krein data of the L2 block of P,
% and the distributed gain is returned as its values on the grid. Write
% P = [Pm Q1; Q2 R] (Pm: R^k -> R^k, Q1: L2 -> R^k, Q2: R^k -> L2,
% R: L2 -> L2), Z = [Z1 Z2] (Z1 nu x k, Z2(s) nu x m), and Rh = R^{-1} in
% the pointwise form (Shivakumar, Das, Peet, arXiv 2208.13104, Lemma 16
% and Cor. 17)
%   Rh0(s) = R0(s)^{-1},
%   Rh(s,t) = M(s)(I - Pg)N(t) (t <= s),   -M(s)Pg N(t) (t > s),
% with M, N and the matrix Pg from the Gohberg-Krein system, solved by RK4
% on a grid of N nodes (GK_GRID). With E(s) = int_a^s Z2 M,
% F(s) = int_a^s Q1 M and Y(s) = int_a^s N Q2, integrated by the same RK4,
%   (Z2 Rh)(s) = Z2(s)Rh0(s) + [(E(b) - E(s))(I - Pg) - E(s)Pg] N(s),
%   (Q1 Rh)(s) = Q1(s)Rh0(s) + [(F(b) - F(s))(I - Pg) - F(s)Pg] N(s),
%   (Rh Q2)(s) = Rh0(s)Q2(s) + M(s)[(I - Pg)Y(s) - Pg(Y(b) - Y(s))],
%   T = Pm - int_a^b Q1 (Rh Q2),   J = int_a^b Z2 (Rh Q2)       (Simpson),
% and the block inverse through the Schur complement T (arXiv 2208.13104,
% Lemma 18) gives the gain
%   u = K1 x1 + int_a^b K2(s) x2(s) ds,
%   K1 = (Z1 - J) T^{-1},    K2(s) = (Z2 Rh)(s) + (J - Z1) T^{-1} (Q1 Rh)(s).
% With no R^k space, K1 is empty and K2 = Z2 Rh. For the 'opvar' form of
% the gain, K2 is fitted on the grid by a polynomial in s (Chebyshev least
% squares), and the relative RMS residual of the fit is returned.
%
% INPUT
% - P:      'copvar' object (or 'opvar', converted by OPVAR2COPVAR), the
%           solved Lyapunov operator: the same spaces and dimensions in and
%           out, one spatial variable, one L2^m[s] space and at most one R^k
%           space; R0(s) invertible at every grid node and T invertible;
% - Z:      'copvar' object (or 'opvar'), the solved free operator
%           R^k x L2^m[s] -> R^nu: one output space R^nu, the input spaces
%           those of P;
% - opts:   (optional) 'struct', each field optional: N (grid nodes,
%           default 101, an even N raised by one), tol_rank (relative SVD
%           cut of the kernel ranks, 1e-13), tol (target relative RMS
%           residual of the fit of K2, 1e-8), deg0 (first fit degree, 4),
%           degmax (largest fit degree, 16).
%
% OUTPUT
% - K:      'struct' with fields
%             P      'double' nu x k, K1;
%             s, w   'double' 1 x N, the grid and its Simpson weights;
%             Q1     'double' nu x m x N, the values K2(s_i);
%             apply  'function_handle', u = apply(x1,x2) = K1 x1 +
%                    sum_i w_i K2(s_i) x2(s_i), x2 an m x N array of values
%                    on K.s or a 'function_handle' s -> m x 1;
%             op     'opvar' object with P = K1 and Q1 the polynomial fit of
%                    K2, R^k x L2^m[s] -> R^nu;
%             c      'copvar' object, K.op converted;
%             fit    'struct' with fields d (the fit degree) and relrms (the
%                    relative RMS residual of the fit);
% - info:   'struct' with fields N, r1 and r2 (kernel ranks), n (= r1 + r2),
%           rcondU22 (reciprocal condition number of the block U22(b) of the
%           Gohberg-Krein solution), condR0 (largest condition number of R0
%           on the grid) and condT (condition number of T; NaN with no R^k
%           space).
%
% See also GETCONTROLLER_SOP, GK_GRID, CHEBFIT_MONO, COPVAR/INV, CLOSEDLOOPPIE.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - getController_direct_sop
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
% Initial coding MMP, 10/09/2026.
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), step comments in the body; code
%                unchanged. The deleted header text was a reference to the
%                delay-case construction it follows, the uses of K.op
%                elsewhere, a comparison with getController_sop with its
%                measurement, and a note on the observer gain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Options: the defaults, overridden by the fields of OPTS
if nargin<3,    opts = struct();    end
o = struct('N',101,'tol',1e-8,'deg0',4,'degmax',16,'tol_rank',1e-13);
fn = fieldnames(opts);
for k = 1:numel(fn),    o.(fn{k}) = opts.(fn{k});   end
% 'opvar' operators to containers; containers in one spatial variable
if isa(P,'opvar'),  P = opvar2copvar(P);    end
if isa(Z,'opvar'),  Z = opvar2copvar(Z);    end
if ~isa(P,'copvar') || ~isa(Z,'copvar')
    error('getController_direct_sop:class','P and Z must be copvar (or opvar) objects.');
end
if numel(P.vars)>1 || numel(Z.vars)>1
    error('getController_direct_sop:nd','One spatial variable only.');
end
% ---- the spaces of P: one L2 space, at most one R^k space, square
[Mp,Np] = size(P);
if Mp~=Np || ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
    error('getController_direct_sop:square','P must map a space list to itself.');
end
l = find(any(P.space_out,2));   r = find(~any(P.space_out,2));
if numel(l)~=1 || numel(r)>1
    error('getController_direct_sop:spaces','P needs one L2[s] space and at most one R^k space.');
end
m = P.dim_out(l);   k = 0;  if ~isempty(r),  k = P.dim_out(r);   end
% ---- the spaces of Z: one output row R^nu, input columns those of P
[Mz,~] = size(Z);
if Mz~=1 || any(Z.space_out(1,:))
    error('getController_direct_sop:Zout','Z must map to R^nu (one finite-dimensional output space).');
end
nu = Z.dim_out(1);
zl = find(any(Z.space_in,2));   zr = find(~any(Z.space_in,2));
if numel(zl)~=1 || Z.dim_in(zl)~=m || (k>0 && (numel(zr)~=1 || Z.dim_in(zr)~=k)) ...
        || (k==0 && ~isempty(zr) && any(Z.dim_in(zr)>0))
    error('getController_direct_sop:Zin','The input spaces of Z must be the output spaces of P (R^%d x L2^%d).',k,m);
end
% ---- block evaluators
% Z1 and Pm as matrices; Z2(s), Q1(s) (left factors) and Q2(s) (right
% factor) as function handles
Z1 = zeros(nu,k);
if k>0 && ~isempty(zr),     Z1 = blockmat(Z.C{1,zr},nu,k);  end
Z2f = eval_RL(Z.C{1,zl},nu,m);
lf = {Z2f};     rf = {};    Pm = zeros(k);
if k>0
    Q1f = eval_RL(P.C{r,l},k,m);    Q2f = eval_LR(P.C{l,r},m,k);    Pm = blockmat(P.C{r,r},k,k);
    lf{2} = Q1f;    rf = {Q2f};
end
% ---- the Gohberg-Krein data and the cumulative integrals
G = gk_grid(P.C{l,l},struct('N',o.N,'tol_rank',o.tol_rank),lf,rf);
N = G.N;    n = G.n;    sg = G.sg;  w = G.w;    Pg = G.Pg;  ImP = eye(n) - Pg;
% (Z2 Rh)(s_i) at the nodes, from E(s) = int_a^s Z2 M
E = G.E{1}; EN = E(:,:,N);
Z2v = zeros(nu,m,N);    K2 = zeros(nu,m,N);
for i = 1:N
    Z2v(:,:,i) = Z2f(sg(i));
    K2(:,:,i) = Z2v(:,:,i)*G.R0inv(:,:,i) + ((EN - E(:,:,i))*ImP - E(:,:,i)*Pg)*G.Nn(:,:,i);   % (Z2 Rh)(s_i)
end
% With an R^k space: (Q1 Rh)(s_i) and (Rh Q2)(s_i) at the nodes, T and J
% by Simpson, then K1 and the term (J - Z1) T^{-1} (Q1 Rh) of K2
K1 = zeros(nu,k);   condT = NaN;
if k>0
    F = G.E{2};     FN = F(:,:,N);  Y = G.Y{1};     YN = Y(:,:,N);
    Q1R = zeros(k,m,N);     T = Pm;     J = zeros(nu,k);
    for i = 1:N
        Q1R(:,:,i) = Q1f(sg(i))*G.R0inv(:,:,i) + ((FN - F(:,:,i))*ImP - F(:,:,i)*Pg)*G.Nn(:,:,i);
        RQ2 = G.R0inv(:,:,i)*Q2f(sg(i)) + G.M(:,:,i)*(ImP*Y(:,:,i) - Pg*(YN - Y(:,:,i)));
        T = T - w(i)*(Q1f(sg(i))*RQ2);
        J = J + w(i)*(Z2v(:,:,i)*RQ2);
    end
    condT = cond(T);
    K1 = (Z1 - J)/T;
    JT = (J - Z1)/T;
    for i = 1:N,    K2(:,:,i) = K2(:,:,i) + JT*Q1R(:,:,i);    end
end
% ---- the one fit: K2(s) as a polynomial for the legacy closed loop
Yk = reshape(permute(K2,[3 1 2]),N,nu*m);
[Cmono,relrms,d] = chebfit_mono(sg,Yk,[G.a G.b],o.deg0,o.degmax,o.tol);
sname = P.vars{1};
Q1poly = polynomial(sparse(Cmono),(0:d)',{sname},[nu m]);
% The 'opvar' gain: P = K1, Q1 = the fitted K2, on the domain and variable
% of P
opvar Kop;
Kop.I = [G.a G.b];
Kop.var1 = polynomial({sname});     Kop.var2 = polynomial({[sname '_dum']});
Kop.P = K1;     Kop.Q1 = Q1poly;
% Outputs
K = struct('P',K1,'s',sg,'w',w,'Q1',K2,'apply',@(x1,x2) apply_gain(K1,K2,sg,w,x1,x2),...
           'op',Kop,'c',opvar2copvar(Kop),'fit',struct('d',d,'relrms',relrms));
info = struct('N',N,'r1',G.r1,'r2',G.r2,'n',n,'rcondU22',G.rcondU22,'condR0',G.condR0,'condT',condT);
end


function u = apply_gain(K1,K2,sg,w,x1,x2)
% U = APPLY_GAIN(K1,K2,SG,W,X1,X2) returns u = K1 x1 + sum_i w_i K2(s_i)
% x2(s_i), the gain applied to the state (x1,x2) by the quadrature (SG,W).
% INPUT  K1:    'double' nu x k;
%        K2:    'double' nu x m x N, the values K2(s_i);
%        sg, w: 'double' 1 x N, the nodes and weights;
%        x1:    'double', k entries;
%        x2:    'double' m x N, the values on SG, or 'function_handle'
%               s -> m x 1.
% OUTPUT u:     'double' nu x 1.
% The values of x2 on the grid
if isa(x2,'function_handle')
    X2 = zeros(size(K2,2),numel(sg));
    for i = 1:numel(sg),    X2(:,i) = x2(sg(i));    end
else
    X2 = x2;
end
% K1 x1 plus the quadrature of K2 x2
u = K1*x1(:);
for i = 1:numel(sg),    u = u + w(i)*(K2(:,:,i)*X2(:,i));    end
end

function Mat = blockmat(B,p,q)
% MAT = BLOCKMAT(B,P,Q) returns the p x q matrix of an R^q -> R^p block,
% an error if its parameter has another size.
% INPUT  B:    'sopvar' object of an R^q -> R^p block, or [] (zero);
%        p, q: 'double', the dimensions.
% OUTPUT Mat:  'double' p x q.
if isempty(B),  Mat = zeros(p,q);   return,     end
Mat = full(B.params{1});
if ~isequal(size(Mat),[p q])
    error('getController_direct_sop:matblock','An R^%d -> R^%d block stores a %d x %d parameter.',q,p,size(Mat,1),size(Mat,2));
end
end

function f = eval_RL(B,p,m)
% F = EVAL_RL(B,P,M) returns the kernel of an L2^m -> R^p block,
% (B x) = int_a^b Bv(t) x(t) dt with Bv(t) = C (I_m kron ZR(t)), as a
% function handle.
% INPUT  B:    'sopvar' object of the block, or [] (zero);
%        p, m: 'double', the dimensions.
% OUTPUT f:    'function_handle', t -> Bv(t), p x m.
if isempty(B),  f = @(t) zeros(p,m);    return,     end
C = B.params{1};    zR = B.ZR{1}(:);
f = @(t) full(C*kron(speye(m),sparse(t.^zR)));
end

function f = eval_LR(B,m,q)
% F = EVAL_LR(B,M,Q) returns the kernel of an R^q -> L2^m block,
% (B v)(s) = Bv(s) v with Bv(s) = (I_m kron ZL(s)') C, as a function
% handle.
% INPUT  B:    'sopvar' object of the block, or [] (zero);
%        m, q: 'double', the dimensions.
% OUTPUT f:    'function_handle', s -> Bv(s), m x q.
if isempty(B),  f = @(s) zeros(m,q);    return,     end
C = B.params{1};    zL = B.ZL{1}(:);
f = @(s) full(kron(speye(m),sparse(s.^zL)')*C);
end
