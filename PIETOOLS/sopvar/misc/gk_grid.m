function G = gk_grid(R,opts,lf,rf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G = GK_GRID(R,OPTS[,LF,RF]) the Gohberg-Krein data of a square 1-D
% 'sopvar' block
%   (R x)(s) = R0(s) x(s) + int_a^s R1(s,t) x(t) dt + int_s^b R2(s,t) x(t) dt
% on a uniform grid of OPTS.N nodes of [a,b], from which R^{-1} is read off
% pointwise (Shivakumar, Das, Peet, arXiv 2208.13104, Lemma 16 and Cor. 17):
%   R^{-1} = {R0^{-1}, M(s)(I - Pg)N(t) for t <= s, -M(s) Pg N(t) for t > s},
%   M(s) = C(s)U(s),  N(t) = V(t)B(t)R0(t)^{-1},
% with C = -R0^{-1}[F1 F2], B = [G1; -G2] from the kernel factors R_i =
% F_i G_i (one SVD of each stored coefficient matrix C_i, the block storing
% R_i = (I kron ZL') C_i (I kron ZR)), U' = BCU, V' = -VBC, U(a) = V(a) = I by
% RK4 with midpoints, and Pg = [0 0; U22(b)^{-1} U21(b) I]. Optionally the
% cumulative integrals a gain construction needs,
%   E_i(s) = int_a^s LF_i(t) M(t) dt,     Y_j(s) = int_a^s N(t) RF_j(t) dt,
% are integrated alongside U and V (one augmented linear RK4 system each),
% so they carry the RK4 accuracy rather than a quadrature of node values.
%
% INPUTS
% - R:     'sopvar', L2^m[s] -> L2^m[s] in one variable, square;
% - opts:  struct, fields N (nodes, default 101, made odd so that Simpson
%          weights exist) and tol_rank (relative SVD cut, 1e-13);
% - lf:    (optional) cell of function handles s -> p_i x m matrix;
% - rf:    (optional) cell of function handles s -> m x q_j matrix.
% OUTPUTS
% - G:     struct: a, b, m, N, sg (1 x N nodes), h, w (1 x N Simpson
%          weights), r1, r2, n = r1 + r2, Pg (n x n), rcondU22, condR0 (max
%          over nodes), R0inv (m x m x N), M (m x n x N), Nn (n x m x N),
%          E (cell, p_i x n x N), Y (cell, n x q_j x N), and the evaluator
%          handles R0fun, Ffun (s -> [F1 F2]), Gfun (s -> [G1; -G2]).
%
% NOTES
% Shared by @sopvar/inv (which fits the pointwise inverse to polynomials)
% and getController_direct_sop (which never forms the inverse). Cost
% O(N n^3) plus O(N n sum p_i) and O(N n sum q_j); nothing scales with q.
% A singular R0 at a node or a near-singular U22(b) is reported through
% condR0 / rcondU22 and a warning.
%
% See also SOPVAR/INV, GETCONTROLLER_DIRECT_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - gk_grid
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
% Initial coding MMP, 10/09/2026: the grid part of @sopvar/inv (same date),
%                moved here so that the gain construction can share it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 || isempty(opts),   opts = struct();    end
if nargin<3,    lf = {};    end
if nargin<4,    rf = {};    end
o = struct('N',101,'tol_rank',1e-13);
fn = fieldnames(opts);
for k = 1:numel(fn),    if isfield(o,fn{k}),    o.(fn{k}) = opts.(fn{k});   end,    end
if mod(o.N,2)==0,   o.N = o.N + 1;  end

% ---- the block must be 1-D, L2 -> L2 in one variable, and square
if numel(R.vars.in)~=1 || numel(R.vars.out)~=1 || ~strcmp(R.vars.in{1},R.vars.out{1})
    error('gk_grid:dim','The block must map L2^m[s] -> L2^m[s] in one variable.');
end
m = R.dims(1);
if R.dims(2)~=m,    error('gk_grid:square','The block must be square (dims %d x %d).',R.dims(1),R.dims(2));    end
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
L = numel(lf);  pL = zeros(1,L);    for i = 1:L,   pL(i) = size(lf{i}(a),1);    end
Rr = numel(rf); qR = zeros(1,Rr);   for j = 1:Rr,  qR(j) = size(rf{j}(a),2);    end
nX = n + sum(pL);   nW = n + sum(qR);       % rows of the left system, columns of the right

% ---- grid and stage data at nodes and midpoints
N = o.N;    sg = linspace(a,b,N);   h = sg(2) - sg(1);    smid = sg(1:end-1) + h/2;
R0inv = zeros(m,m,N);   condR0 = zeros(N,1);
Cn = cell(N,1);     Bn = cell(N,1);
Lnode = cell(N,1);  Lmid = cell(N-1,1);     % left system:  X' = Lmat * X(1:n,:), X = [U; E_1; ...]
Rnode = cell(N,1);  Rmid = cell(N-1,1);     % right system: W' = W(:,1:n) * Rmat, W = [V, Y_1, ...]
for i = 1:N
    [Li,Ri] = bases(sg(i));
    R0i = full(Li*C0*Ri);   condR0(i) = cond(R0i);
    if rcond(R0i)<eps
        % a PI operator whose multiplier is zero or singular at a point has no
        % bounded inverse on L2 (compact, or 0 in the essential spectrum); the
        % Gohberg-Krein formulas need R0(s)^{-1} at every node. For the QT
        % synthesis form use the xi-coordinates (map Sec. 8.2) instead.
        error('gk_grid:multiplier',['The multiplier R0(s) is singular at s = %.4g (rcond %.2g): the operator ' ...
              'is not boundedly invertible on L2 and the Gohberg-Krein construction does not apply.'],sg(i),rcond(R0i));
    end
    R0inv(:,:,i) = R0i\eye(m);
    [Cn{i},Bn{i}] = CB(Li,Ri,R0inv(:,:,i));
    [Lnode{i},Rnode{i}] = stage(sg(i),Cn{i},Bn{i},R0inv(:,:,i));
end
for i = 1:N-1
    [Li,Ri] = bases(smid(i));
    R0i = full(Li*C0*Ri);   R0inv_i = R0i\eye(m);
    [Ci,Bi] = CB(Li,Ri,R0inv_i);
    [Lmid{i},Rmid{i}] = stage(smid(i),Ci,Bi,R0inv_i);
end

% ---- RK4 on the two augmented linear systems
X = zeros(nX,n,N);  W = zeros(n,nW,N);
X(1:n,:,1) = eye(n);    W(:,1:n,1) = eye(n);
for k = 1:N-1
    Xk = X(:,:,k);  Wk = W(:,:,k);
    k1 = Lnode{k}*Xk(1:n,:);    k2 = Lmid{k}*(Xk(1:n,:) + h/2*k1(1:n,:));
    k3 = Lmid{k}*(Xk(1:n,:) + h/2*k2(1:n,:));   k4 = Lnode{k+1}*(Xk(1:n,:) + h*k3(1:n,:));
    X(:,:,k+1) = Xk + h/6*(k1 + 2*k2 + 2*k3 + k4);
    k1 = Wk(:,1:n)*Rnode{k};    k2 = (Wk(:,1:n) + h/2*k1(:,1:n))*Rmid{k};
    k3 = (Wk(:,1:n) + h/2*k2(:,1:n))*Rmid{k};   k4 = (Wk(:,1:n) + h*k3(:,1:n))*Rnode{k+1};
    W(:,:,k+1) = Wk + h/6*(k1 + 2*k2 + 2*k3 + k4);
end
U = X(1:n,:,:);     V = W(:,1:n,:);

% ---- Pg from U(b): the upper kernel's block decides invertibility
Ub = U(:,:,N);
U21 = Ub(r1+1:n,1:r1);  U22 = Ub(r1+1:n,r1+1:n);
if r2>0,    rc = rcond(U22);    else,   rc = 1;     end
if rc<1e-14
    warning('gk_grid:singular','U22(b) is near singular (rcond %.2g): the block is likely not invertible.',rc);
end
Pg = zeros(n);
if r2>0,    Pg(r1+1:n,1:r1) = U22\U21;  Pg(r1+1:n,r1+1:n) = eye(r2);    end

% ---- node values of M, N and the cumulative integrals
M = zeros(m,n,N);   Nn = zeros(n,m,N);
for i = 1:N
    M(:,:,i) = Cn{i}*U(:,:,i);
    Nn(:,:,i) = V(:,:,i)*Bn{i}*R0inv(:,:,i);
end
E = cell(1,L);  off = n;
for i = 1:L,    E{i} = X(off+(1:pL(i)),:,:);    off = off + pL(i);  end
Y = cell(1,Rr); off = n;
for j = 1:Rr,   Y{j} = W(:,off+(1:qR(j)),:);    off = off + qR(j);  end
w = (h/3)*[1, repmat([4 2],1,(N-3)/2), 4, 1];                              % composite Simpson, N odd

G = struct('a',a,'b',b,'m',m,'N',N,'sg',sg,'h',h,'w',w,'r1',r1,'r2',r2,'n',n,...
           'Pg',Pg,'rcondU22',rc,'condR0',max(condR0),'R0inv',R0inv,'M',M,'Nn',Nn,...
           'R0fun',@(s) R0_at(s),'Ffun',@(s) F_at(s),'Gfun',@(s) G_at(s));
G.E = E;    G.Y = Y;        % assigned after: an empty cell inside struct() makes an empty struct array

    function [Lb,Rb] = bases(s)
        Lb = kron(Im,sparse(s.^zL)');   Rb = kron(Im,sparse(s.^zR));
    end
    function [Ci,Bi] = CB(Lb,Rb,R0i)
        % C(s) = -R0(s)^{-1} [F1 F2], B(s) = [G1; -G2]
        Fm = [full(Lb*U1), full(Lb*U2)];
        Gm = [full(V1*Rb); -full(V2*Rb)];
        Ci = -R0i*Fm;   Bi = Gm;
    end
    function [Lmat,Rmat] = stage(s,Ci,Bi,R0i)
        A = Bi*Ci;
        Lmat = zeros(nX,n);     Lmat(1:n,:) = A;    off1 = n;
        for ii = 1:L
            Lmat(off1+(1:pL(ii)),:) = lf{ii}(s)*Ci;     off1 = off1 + pL(ii);
        end
        Rmat = zeros(n,nW);     Rmat(:,1:n) = -A;   off2 = n;
        for jj = 1:Rr
            Rmat(:,off2+(1:qR(jj))) = Bi*R0i*rf{jj}(s);     off2 = off2 + qR(jj);
        end
    end
    function R0s = R0_at(s)
        [Lb,Rb] = bases(s);     R0s = full(Lb*C0*Rb);
    end
    function Fs = F_at(s)
        [Lb,~] = bases(s);      Fs = [full(Lb*U1), full(Lb*U2)];
    end
    function Gs = G_at(s)
        [~,Rb] = bases(s);      Gs = [full(V1*Rb); -full(V2*Rb)];
    end
end


function C = zero_or(C,nr,nc)
% [] and scalar 0 are the zero-block shorthand of the class
if isempty(C) || (isscalar(C) && C==0),    C = sparse(nr,nc);   end
if ~isequal(size(C),[nr nc])
    error('gk_grid:params','A parameter is %d x %d; the bases give %d x %d.',size(C,1),size(C,2),nr,nc);
end
end

function [U,V] = factor_C(C,tol)
% C = U V with the numerical rank: economy SVD, singular values below tol*max dropped
if nnz(C)==0,   U = zeros(size(C,1),0);  V = zeros(0,size(C,2));  return,  end
[Us,S,Vs] = svd(full(C),'econ');
sv = diag(S);   r = find(sv>=tol*sv(1),1,'last');
U = Us(:,1:r)*diag(sv(1:r));    V = Vs(:,1:r)';
end
