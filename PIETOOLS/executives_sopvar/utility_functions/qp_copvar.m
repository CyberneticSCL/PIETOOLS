function [X,info] = qp_copvar(A,B,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,INFO] = QP_COPVAR(A,B[,OPTS]) solves the least squares problem
%
%   min_X || A X - B ||^2_L2     over X with every block a polynomial of
%                                degree <= deg in each spatial variable,
%
% for a 'copvar' A (M x N grid of 'sopvar' blocks) and a 'copvar' B (M x 1,
% blocks R -> L2, no input variables), as a quadratic program in the
% coefficients of X, which are the decision variables.
%
% Construction. Every block j of X is linear in its coefficient vector x_j,
%   X_j = Xh_j x_j,    Xh_j: R^{n_j} -> L2^{p_j}[s^j]  with params = I,
% a 'sopvar' on a monomial basis with degrees 0..deg in each variable
% (n_j = p_j*prod(deg+1), n_j = p_j on an R^p space). With the container Xh
% = blkdiag(Xh_j) the objective is
%   ||A X - B||^2 = x' H x - 2 f' x + c,
%   H = (A Xh)'(A Xh),   f = (A Xh)'B,   c = B'B,
% where A*Xh, (A Xh)'(A Xh), (A Xh)'B, B'B are formed by the container
% composition and the R^n -> R^n, R -> R^n, R -> R blocks of H, f, c are read
% off their parameters (H = Xh'(A'A)Xh as before, without forming A'A). The degree bound is the basis: a coefficient outside it is
% not a variable, which is the constraint "degree of X <= deg".
%
% INPUTS
% - A:     'copvar', M x N;
% - B:     'copvar', M x 1, no input variables (blocks R -> L2^{q_i}[s^i]);
% - opts:  (optional) struct with any of the fields
%   .deg      maximal degree in each spatial variable of X: a scalar (all
%             variables) or a vector with one entry per variable of A.vars.
%             Default 4;
%   .basis    'monomial' (default), 'legendre' or 'chebyshev': the polynomials
%             whose coefficients are the unknowns, per variable and
%             degree 0..deg on the domain of the variable (Legendre
%             orthonormal in L2, Chebyshev of the first kind T_k). The
%             function X is the same set of polynomials in every case; only
%             the conditioning of H changes. H is still assembled through the
%             monomial arithmetic of sopvar, so a rounding error of about
%             1e-16*norm(S)^2 enters (S: monomial coefficients of the basis),
%             which grows quickly with the degree;
%   .assembly how H, f, c are formed: 'container' (default) from the
%             container products AXh = A*Xh, H = AXh'*AXh, f = AXh'*B, B'*B,
%             which carry the monomial coefficients of the basis through
%             sopvar mtimes; or 'apply' from A*phi_u, computed by apply_sopvar
%             for every basis function phi_u and integrated by Gauss-Legendre
%             quadrature, H being the Gram matrix of these functions
%             (positive semidefinite by construction). 'container' is faster
%             and, since A'*A is not formed, light on memory; 'apply' is slower
%             (one apply_sopvar per unknown and block) and needs L2 input
%             spaces, but its H stays accurate at high degree: through the
%             monomial coefficients of the basis the 'container' H loses
%             accuracy from degree 6 in 2-D (indefinite from degree 6 to 8,
%             1e4 to 1e11 negative eigenvalues at degrees 8 to 10), from about
%             degree 10 to 12 in 1-D. Measured, 2-D, 2 x 2 blocks, degrees 2 to
%             10 (36 to 484 unknowns): 'container' 0.7 to 3.5 s and 0.9 to 1.6 GB
%             of peak memory (1.9 to 121 s and 1.1 to 36 GB when H was formed as
%             Xh'*(A'*A)*Xh), 'apply' 11 to 309 s and 0.95 to 1.7 GB. Use
%             'container' up to degree 4 to 6 and 'apply' above;
%   .solver   'quadprog' (the one on the path: MATLAB's, or the MOSEK
%             toolbox's, whose options are set with mskoptimset; the MOSEK
%             case is untested), 'backslash' (the normal equations through
%             LSQMINNORM) or 'auto' (MATLAB's quadprog if that is the one on
%             the path, else backslash). The QP has no constraints, and
%             backslash was the more accurate of the two on the problems
%             tested (quadprog: 1e-6 to 1e-4 against 1e-10 for the error in
%             X on ill-conditioned H). Default 'backslash';
%   .lambda   Tikhonov term lambda*||x||^2 added to the objective, >= 0.
%             Default 0;
%   .verbose  print the sizes and the residual. Default false.
%   Unknown fields are an error.
%
% RECOMMENDATION (measured, 2-D, 2 x 2 blocks, deg 2 to 10; 1-D similar)
% - The residual ||A X - B||^2/||B||^2 of the returned X is about the same for
%   all bases up to degree 6, and at degree 10 it is at the 1e-9 level in all
%   three with 'container'. The choices differ in H and in the accuracy of x:
% - 'monomial': H is nearly singular already at degree 2 (cond 1e9, 1e19 to
%   1e21 from degree 6), whatever the assembly. Use only for deg <= 4.
% - 'chebyshev': no advantage over 'legendre'; its H is worse conditioned in L2
%   (cond 1e4 to 1e5) and it breaks down in the same way with 'container'.
% - 'legendre': H is well conditioned (cond about 5e3 in 2-D) if assembled
%   without monomial coefficients, i.e. with 'apply'.
% - So, deg <= 4 (6 in 2-D): 'legendre' with 'container' (fast, light, same
%   residual); higher degree, L2 input spaces: 'legendre' with 'apply'
%   (H positive definite, most accurate, but 10 to 300 s in 2-D); take the
%   residual from info.resq and the coefficients from info.x, since the
%   monomial coefficients of the returned X lose accuracy beyond degree 8 to 10.
%   For an R^p input space only 'container' is available, so keep deg low.
%
% OUTPUTS
% - X:     'copvar', N x 1, blocks R -> L2^{p_j}[s^j] on the monomial basis
%          of degrees 0..deg, the minimizer;
% - info:  struct with fields
%   .n        number of coefficients (decision variables),
%   .obj      x'Hx - 2f'x + c = ||A X - B||^2 (cancellation limits it to
%             about 1e-16*c, so a relative residual below 1e-8 is noise),
%   .relres   sqrt(max(obj,0)/c) = ||A X - B||/||B||,
%   .resq     ||A X - B||/||B|| from the quadrature samples, assembly 'apply'
%             only (NaN otherwise): a sum of squares of samples, so unlike
%             .relres it has no cancellation floor, and it is the residual of
%             the coefficients x, not of the stored monomial form of X;
%   .x        the minimizer's coefficients in the chosen basis (the unknowns).
%             For large degrees they can be huge (1e13 and more) and X, which
%             stores them as monomial coefficients, is then numerically useless
%             even when x is right: judge the fit by .resq;
%   .solver   the solver used,  .exitflag  its exit flag,
%   .cond     condition number of H (n <= 2000 only, else NaN).
%
% NOTES
% H is the Gram matrix of the monomials under A'A, Hilbert-like, so it is
% numerically singular for more than about 10 degrees per variable. The
% unknowns are scaled by sqrt(diag(H)) before the solve, and a nonzero lambda
% or the 'backslash' minimum norm solution regularizes the rest. The number
% of unknowns grows as prod(deg+1) per block, so the dense H limits N-D use.
%
% See also COPVAR/MTIMES, SOPVAR.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - qp_copvar
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
% Initial coding AT - 10/26/2026 

if nargin<3 || isempty(opts),  opts = struct();  end
%o = struct('deg',4,'basis','monomial','assembly','container','solver','auto','lambda',0,'verbose',false);   % AT, 10/28/2026 (was)
o = struct('deg',4,'basis','monomial','assembly','container','solver','backslash','lambda',0,'verbose',false);   % AT, 10/28/2026
for f = fieldnames(opts)'
    if ~isfield(o,f{1}),  error('qp_copvar:opts','Unknown option ''%s''.',f{1});  end
    if ~isempty(opts.(f{1})),  o.(f{1}) = opts.(f{1});  end
end
if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('A and B should be copvar objects.')
end
if size(B,2)~=1 || size(B,1)~=size(A,1)
    error('B must be an M x 1 copvar with M = size(A,1).')
end
if ~isscalar(o.lambda) || o.lambda<0,  error('lambda must be a scalar >= 0.');  end
solver = lower(o.solver);
if strcmp(solver,'auto')
    % The quadprog on the path may be the one of the MOSEK toolbox, which takes   % AT, 10/28/2026
    % 'mskoptimset' options, not 'optimoptions': detect it by its location, only  % AT, 10/28/2026
    % for the solver 'auto' (the solver 'quadprog' does it where it calls it).    % AT, 10/28/2026
    qp_file = which('quadprog');                                            % AT, 10/28/2026
    qp_is_mosek = ~isempty(qp_file) && contains(lower(qp_file),'mosek');    % AT, 10/28/2026
%    if exist('quadprog','file'),  solver = 'quadprog';  else,  solver = 'backslash';  end   % AT, 10/28/2026 (was)
    % 'quadprog' of the Optimization Toolbox if it is the one on the path, else backslash  % AT, 10/28/2026
    if ~isempty(qp_file) && ~qp_is_mosek,  solver = 'quadprog';  else,  solver = 'backslash';  end   % AT, 10/28/2026
end
if ~any(strcmp(solver,{'quadprog','backslash'})),  error('qp_copvar:opts','Unknown solver ''%s''.',solver);  end
assembly = lower(o.assembly);
if ~any(strcmp(assembly,{'container','apply'})),  error('qp_copvar:opts','Unknown assembly ''%s''.',assembly);  end
basis = lower(o.basis);
if ~any(strcmp(basis,{'monomial','legendre','chebyshev'})),  error('qp_copvar:opts','Unknown basis ''%s''.',basis);  end
N = size(A,2);
nv = numel(A.vars);
deg = o.deg(:)';
if isscalar(deg),  deg = repmat(deg,1,nv);  end
if numel(deg)~=nv || any(deg<0) || any(deg~=round(deg))
    error('qp_copvar:deg','deg must be a nonnegative integer scalar, or one per variable of A.vars (%d).',nv);
end

% QP in the coefficients x of X: min x'Hx - 2f'x + c  (= ||A X - B||^2):
%   N                            number of blocks of X (columns of A)
%   nv                           number of spatial variables in A.vars
%   deg                          1 x nv maximal degree of X in each variable
%   mask                         1 x nv, the variables of block j of X (A.space_in(j,:))
%   p                            components of block j (A.dim_in(j))
%   vnames, dj                   names and maximal degrees of those variables
%   ZL                           monomial exponents 0..dj(k) of block j, per variable
%   nZ                           monomials per component, prod(dj+1) (1 on an R^p space)
%   nj                           N x 1 coefficients of block j, p*nZ
%   off, n                       offsets of the blocks in x, n = sum(nj) unknowns
%   S_k, Sj                      monomial coefficients of the basis functions of
%                                variable k, and kron(S_1,...,S_nv) for block j
%   Xh                           container, Xh_j: R^{nj} -> L2^p[s^j], params kron(I_p,Sj),
%                                so X_j = Xh_j x_j is linear in x_j
%   AXh                          A*Xh, container, blocks R^{nj} -> L2^{q_i}[s^i]
%   Hc, fc, cc                   AXh'*AXh (R^nj -> R^ni blocks), AXh'*B (R -> R^ni),
%                                B'*B (R -> R), as containers
%   H, f, c                      n x n matrix, n x 1 vector and scalar read off the
%                                parameters of Hc, fc, cc
%   s, Hy, fy, y                 scaling sqrt(diag(H)), scaled H and f, scaled unknowns
%                                y = x.*s (conditioning only)
%   x                            the minimizer, the coefficients of X
%   obj                          x'Hx - 2f'x + c, the squared residual ||A X - B||^2
%   Cx                           N x 1 cell of the blocks R -> L2 of X, params Sj*x_j

% ---- the linear parametrization: block j of X is Xh_j x_j, params = I
Ch = cell(N,N);   nj = zeros(N,1);
for j = 1:N
    mask = A.space_in(j,:);   p = A.dim_in(j);
    vnames = A.vars(mask);    dj = deg(mask);
    ZL = cell(1,numel(vnames));
    for k = 1:numel(vnames),  ZL{k} = (0:dj(k))';  end
    nZ = prod(dj+1);                                % 1 on an R^p space
    nj(j) = p*nZ;
    vs = struct('in',{{}},'out',{vnames(:)'});
    ds = struct('in',zeros(0,2),'out',A.dom(mask,:));
    % params = monomial coefficients of the basis functions: kron(I_p,S_1,...,S_nv),
    % first variable outermost as in the monomial index; = I for 'monomial'.
    Sj = 1;   vidx = find(mask);
    for k = 1:numel(vnames)
        if strcmp(basis,'monomial'),  Sk = speye(dj(k)+1);
        else,  Sk = basis_coef(dj(k),A.dom(vidx(k),:),basis);
        end
        Sj = kron(Sj,Sk);
    end
    Ch{j,j} = sopvar({kron(speye(p),sparse(Sj))},vs,ZL,{},ds,[p,nj(j)]);
end
Xh = copvar(Ch);

% ---- H = Xh'(A'A)Xh, f = Xh'(A'B), c = B'B
off = [0;cumsum(nj)];   n = off(end);
if strcmp(assembly,'apply')
    [H,f,c,Qs] = assemble_apply(A,B,Xh,nj,off);     % Gram matrix of A*phi_u by quadrature
else
%    AtA = A'*A;   AtB = A'*B;                                              % AT, 10/28/2026 (was)
%    Hc = Xh'*AtA*Xh;   fc = Xh'*AtB;   cc = B'*B;                          % AT, 10/28/2026 (was)
    % H = (A Xh)'(A Xh), f = (A Xh)'B: A*Xh first, so that A'*A, with twice the    % AT, 10/28/2026
    % monomials of A, is never formed (92% of the allocation of the old product).  % AT, 10/28/2026
    AXh = A*Xh;                                                             % AT, 10/28/2026
    Hc = AXh'*AXh;   fc = AXh'*B;   cc = B'*B;                              % AT, 10/28/2026
    H = zeros(n);   f = zeros(n,1);
    for i = 1:N
        ri = off(i)+1:off(i+1);
        for j = 1:N
            blk = Hc.C{i,j};
            if ~isempty(blk),  H(ri,off(j)+1:off(j+1)) = full(blk.params{1});  end
        end
        blk = fc.C{i,1};
        if ~isempty(blk),  f(ri) = full(blk.params{1});  end
    end
    c = 0;   blk = cc.C{1,1};
    if ~isempty(blk),  c = full(blk.params{1});  end
end
H = (H+H')/2;
Hraw = H;
if o.lambda>0,  H = H + o.lambda*speye(n);  end

% ---- the QP, in the variables y = x.*s with s = sqrt(diag(H))
s = sqrt(diag(H));   s(s==0) = 1;
Hy = H./(s*s');   fy = f./s;   Hy = (Hy+Hy')/2;
exitflag = 0;
switch solver
    case 'quadprog'
        if ~exist('qp_is_mosek','var')      % not set by the solver 'auto': look at the quadprog on the path   % AT, 10/28/2026
            qp_is_mosek = contains(lower(which('quadprog')),'mosek');       
        end                                                                 
        if qp_is_mosek
            qo = mskoptimset('Display','off');   % MOSEK's quadprog: its own options, default tolerances  % AT, 10/28/2026
        else                                                                
            qo = optimoptions('quadprog','Display','off','OptimalityTolerance',1e-12,'StepTolerance',1e-14,'ConstraintTolerance',1e-12);
        end                                                                 
        [y,~,exitflag] = quadprog(2*Hy,-2*fy,[],[],[],[],[],[],[],qo);
    case 'backslash'
        y = lsqminnorm(Hy,fy); % Hy \ fy
end
x = y./s;
obj = x'*H*x - 2*f'*x + c;
info = struct('n',n,'obj',obj,'relres',sqrt(max(obj,0)/max(c,eps)),'solver',solver,'exitflag',exitflag,'cond',NaN,'rank',NaN,'mineig',NaN,'nneg',NaN,'H',[],'f',[]);
if n<=2000,  info.H = Hraw;   info.f = f;  end      % Hraw: H before the lambda term
info.x = x;                                         % the coefficients in the chosen basis (blocks at off, sizes nj)
info.resq = NaN;                                    % ||A X - B||/||B|| from the quadrature samples (assembly 'apply'):
if strcmp(assembly,'apply')                         % sum_i |W_i x - b_i|_w^2, no cancellation in x'Hx - 2f'x + c
    rr = 0;
    for i = 1:numel(Qs.W),  e = Qs.W{i}*x - Qs.b{i};   rr = rr + e'*(Qs.w{i}.*e);  end
    info.resq = sqrt(rr/max(c,eps));
end
if n<=2000                              % cond and numerical rank of H (singular values > n*eps*max)
    sv = svd(H);   info.cond = sv(1)/sv(end);   info.rank = nnz(sv>n*eps(sv(1)));
    ev = eig(H);   info.mineig = min(ev)/max(ev);   info.nneg = nnz(ev<-n*eps(max(ev)));   % H is PSD in exact arithmetic
end
if o.verbose
    fprintf('qp_copvar: %d coefficients, solver %s, ||AX-B||/||B|| = %.3e\n',n,solver,info.relres);
end

% ---- X from the coefficients: the same bases, params = x_j
Cx = cell(N,1);
for j = 1:N
    b = Xh.C{j,j};
    Cx{j} = sopvar({sparse(b.params{1}*x(off(j)+1:off(j+1)))},b.vars,b.ZL,{},b.dom,[b.dims(1),1]);
end
X = copvar(Cx);
end

function [H,f,c,Q] = assemble_apply(A,B,Xh,nj,off)
% H(u,v) = <A phi_u, A phi_v>, f(u) = <A phi_u, B>, c = <B,B> for the basis
% functions phi_u of X (the columns of Xh), with every A phi_u computed by
% apply_sopvar and the L2 inner products by Gauss-Legendre quadrature on the
% output space of each block row. H = W'*diag(weights)*W is a Gram matrix of
% sampled functions, so it is positive semidefinite whatever the basis, and no
% product of the monomial coefficients of the basis enters.
M = size(A,1);   N = size(A,2);   n = off(end);   q = A.dim_out;
Wp = cell(M,n);                                     % Wp{i,u} = polynomial q_i x 1, A_i phi_u
for j = 1:N
    b = Xh.C{j,j};   p = b.dims(1);   nZ = nj(j)/p;
    vn = b.vars.out;
    if isempty(vn),  error('qp_copvar:apply','assembly ''apply'' needs L2 input spaces (an R^p space has no variables).');  end
    Sj = full(b.params{1}(1:nZ,1:nZ));              % kron(S_1,...,S_nv), basis coefficients
    dm = zeros(1,0);   nprev = 1;                   % exponent tuples of the monomial index
    for k = 1:numel(vn)
        ek = b.ZL{k}(:);   dm = [repelem(dm,numel(ek),1),repmat(ek,nprev,1)];   nprev = nprev*numel(ek);
    end
    for cc = 1:p
        for k = 1:nZ
            u = off(j)+(cc-1)*nZ+k;
            phi = polynomial(sparse((1:nZ)',cc*ones(nZ,1),Sj(:,k),nZ,p),dm,vn,[p 1]);   % nterms x entries
            for i = 1:M
                if ~isempty(A.C{i,j}),  Wp{i,u} = apply_sopvar(A.C{i,j},phi);  end
            end
        end
    end
end
H = zeros(n);   f = zeros(n,1);   c = 0;
Q = struct('W',{cell(M,1)},'b',{cell(M,1)},'w',{cell(M,1)});   % samples, for the residual of a solution
for i = 1:M
    mask = A.space_out(i,:);   vr = A.vars(mask);   dom = A.dom(mask,:);   nv = numel(vr);
    Bi = [];   if ~isempty(B.C{i,1}),  Bi = apply_sopvar(B.C{i,1},polynomial(1));  end
    % quadrature exact for the products: nq = maxdeg+1 nodes per variable
    md = zeros(1,nv);
    allp = [Wp(i,:),{Bi}];
    for t = 1:numel(allp)
        w = allp{t};   if isempty(w) || ~isa(w,'polynomial'),  continue,  end
        for k = 1:numel(w.varname)
            col = find(strcmp(vr,w.varname{k}),1);
            if ~isempty(col) && ~isempty(w.degmat),  md(col) = max(md(col),full(max(w.degmat(:,k))));  end
        end
    end
    nq = md+1;   xs = cell(1,nv);   ws = cell(1,nv);
    for k = 1:nv
        [xk,wk] = gauss_legendre(nq(k));   xs{k} = dom(k,1)+(dom(k,2)-dom(k,1))*(xk+1)/2;   ws{k} = wk*(dom(k,2)-dom(k,1))/2;
    end
    [G{1:max(nv,1)}] = ndgrid(xs{:});   [Gw{1:max(nv,1)}] = ndgrid(ws{:});
    pts = zeros(numel(G{1}),nv);   wq = ones(numel(G{1}),1);
    for k = 1:nv,  pts(:,k) = G{k}(:);   wq = wq.*Gw{k}(:);  end
    npt = size(pts,1);   Wr = zeros(npt*q(i),n);
    for u = 1:n
        if ~isempty(Wp{i,u}),  Wr(:,u) = reshape(poly_vals(Wp{i,u},vr,pts,q(i)),[],1);  end
    end
    wts = repmat(wq,q(i),1);
    H = H + Wr'*(wts.*Wr);
    bv = zeros(npt*q(i),1);
    if ~isempty(Bi)
        bv = reshape(poly_vals(Bi,vr,pts,q(i)),[],1);
        f = f + Wr'*(wts.*bv);   c = c + bv'*(wts.*bv);
    end
    Q.W{i} = Wr;   Q.b{i} = bv;   Q.w{i} = wts;
end
end

function V = poly_vals(w,vr,pts,nr)
% values of the polynomial w (nr x 1) at the points pts (columns = variables vr): npt x nr
if ~isa(w,'polynomial'),  V = double(w).*ones(size(pts,1),nr);  return,  end
dmat = w.degmat;   cf = w.coefficient;
Mon = ones(size(pts,1),size(dmat,1));
for k = 1:numel(w.varname)
    col = find(strcmp(vr,w.varname{k}),1);
    if isempty(col),  error('qp_copvar:apply','Variable %s of a result is not an output variable.',w.varname{k});  end
    Mon = Mon.*(pts(:,col).^(full(dmat(:,k)).'));
end
V = Mon*full(cf);
end

function [x,w] = gauss_legendre(n)
% nodes and weights on [-1,1] (Golub-Welsch)
if n==1,  x = 0;  w = 2;  return,  end
k = 1:n-1;   bb = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(bb,1)+diag(bb,-1));
[x,io] = sort(diag(D));   w = 2*V(1,io)'.^2;
end

function S = basis_coef(d,dom,type)
% S(e+1,k+1) = coefficient of s^e in the k-th basis polynomial, k = 0..d, on
% dom = [a,b]: 'legendre' orthonormal in L2(a,b), 'chebyshev' the T_k of the
% first kind in t = (2s-a-b)/(b-a). Built from the three-term recurrence in
% t, then composed with t = alpha*s + beta, in ascending powers.
a = dom(1);   b = dom(2);   al = 2/(b-a);   be = -(a+b)/(b-a);
Pt = zeros(d+1);                                    % column k: coefficients of p_k in t
Pt(1,1) = 1;
if d>=1,  Pt(2,2) = 1;  end
for k = 1:d-1
    shifted = [0; Pt(1:d,k+1)];                     % t*p_k
    if strcmp(type,'legendre')
        Pt(:,k+2) = ((2*k+1)*shifted - k*Pt(:,k))/(k+1);
    else
        Pt(:,k+2) = 2*shifted - Pt(:,k);
    end
end
S = zeros(d+1);
for k = 0:d
    q = zeros(d+1,1);   base = 1;                   % base = (alpha*s+beta)^i
    for i = 0:k
        q(1:numel(base)) = q(1:numel(base)) + Pt(i+1,k+1)*base(:);
        base = conv(base,[be al]);
    end
    S(:,k+1) = q;
end
if strcmp(type,'legendre'),  S = S.*sqrt((2*(0:d)+1)/(b-a));  end
end
