function [Pinv,info] = inv(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [PINV,INFO] = INV(P[,OPTS]) the inverse of a square 'copvar'.
%
% Two or more spatial variables: P must be separable, P = M + F with M a
% block-diagonal multiplier and F a sum of full integral operators (no
% multiplier between different spaces, no term that multiplies in one
% shared variable and integrates in another, and the lower and upper
% kernels of each shared variable equal). Each block of F is the kernel
% (I x Z_L(s))' C (I x Z_R(t)) integrated over the whole domain, so
% F = L H R with L the monomial multipliers, H = blkdiag(C) and R the
% integrals of the monomials, and by the Woodbury identity
%   P^{-1} = M^{-1} + M^{-1} L Hh R M^{-1},   Hh = -H (I + K H)^{-1},
%   K = R M^{-1} L,
% K by Gauss-Legendre quadrature with M^{-1} inverted at each node, and the
% factors M^{-1} [I, L] and R M^{-1} fitted by tensor Chebyshev least squares
% at degree d per variable, raised from deg0 by 2 until the relative RMS
% residual is below tol or d = degmax.
%
% One spatial variable: P over the spaces R^k and L2^m[s] (either may be
% absent), the container counterpart of @opvar/inv. With the grid written as
%   P = [Pm  Q1;  Q2  R],   Pm: R^k -> R^k,  Q1: L2 -> R^k,  Q2: R^k -> L2,
% R: L2 -> L2, the L2 block is inverted by @sopvar/inv (Gohberg-Krein) and
% the rest is the block inverse through the finite-dimensional Schur
% complement (Shivakumar, Das, Peet, arXiv 2208.13104, Lemma 18):
%   T = Pm - Q1 Rh Q2  (k x k),   Rh = R^{-1},
%   P^{-1} = [T^{-1},        -T^{-1} Q1 Rh;
%             -Rh Q2 T^{-1},  Rh + Rh Q2 T^{-1} Q1 Rh].
% Only one operator inverse (Rh) and one matrix inverse are taken; the stock
% 'inv_opvar' takes the other Schur complement, (R - Q2 Pm^{-1} Q1)^{-1},
% which is a second operator inverse of a modified kernel. The products are
% block compositions of the class ('sopvar' mtimes), exact on polynomials.
%
% INPUTS
% - P:      'copvar', square: the same spaces and dimensions in and out; in
%           one spatial variable at most one R^k and one L2[s] space, in
%           several any spaces, P separable;
% - opts:   (optional) 'struct': in one variable the options of @sopvar/inv
%           (N, tol, deg0, degmax, deg, tol_rank); in several, tol (1e-8),
%           deg0 (4), degmax (14), deg (one fixed degree), nq (Gauss nodes
%           per variable for K, 40) and septol (relative tolerance of the
%           separability test, 1e-10).
% OUTPUTS
% - Pinv:   'copvar' over the same spaces in the same order as P;
% - info:   one variable: the INFO of @sopvar/inv (empty when there is no
%           L2 space), plus condT, the condition number of the Schur
%           complement, when both spaces are present. Several variables:
%           'struct' with fields method ('separable'), d (fit degree),
%           relrms (largest relative RMS residual of the fits), condIKH
%           (condition number of I + K H), Kdiff (relative change of K from
%           nq to nq + 10 nodes), nF (number of integral blocks) and note.
%
% NOTES
% A structurally zero block ([]) of P is respected: a zero Q1 or Q2 leaves
% T = Pm and the off-diagonal blocks of the inverse zero. The R^k -> R^k
% block Pm may itself be [], in which case T = -Q1 Rh Q2.
% In several variables a non-separable P raises an error (no inverse
% formula for a non-separable kernel exists in the toolbox; 'opvar2d' falls
% back to least squares, 'mrdivide', and GETOBSERVER_SOP has the method
% 'qp'). M^{-1} is rational, so the inverse is a polynomial approximation;
% info.relrms and info.d say how good.
%
% See also SOPVAR/INV, GETCONTROLLER_SOP, GETOBSERVER_SOP, COPVAR2OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - inv (copvar)
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
% MMP, 11/04/2026: Two or more spatial variables: the inverse of a separable
%                P by the Woodbury identity (local INV_SEP), the container
%                counterpart of the stock inv_opvar2d_separable, in any
%                number of variables and spaces; before, N-D containers
%                were refused. The 2-D H-infinity estimator gain uses it in
%                place of the stock getObserver_2D. Header: the N-D method,
%                its options and info fields.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): INPUT/OUTPUT for
%                the local functions 'blockmat' and 'matblock'; code
%                unchanged.
% MMP, 11/09/2026: Inlining pass (maintainer, 11/09/2026: locals that are
%                very short and called once go into the body): cheb1 into
%                cheb_fit, cheb_grid into fits, gl_grid (with gauss1 inside
%                it) into kmat, each with its header as the comment of the
%                block; tensor_points stays, called from both grids.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2,    opts = struct();    end
if ~isa(P,'copvar'),    error('copvar:inv:class','Input must be a copvar.');   end
% if numel(P.vars)>1                                                        % MMP, 11/04/2026 (was)
%     error('copvar:inv:nd','INV is implemented for one spatial variable; this container has %d (%s).',numel(P.vars),strjoin(P.vars,',')); % MMP, 11/04/2026 (was)
% end                                                                       % MMP, 11/04/2026 (was)
[M,Nc] = size(P);
if M~=Nc || ~isequal(P.space_out,P.space_in) || ~isequal(P.dim_out(:),P.dim_in(:))
    error('copvar:inv:square','The container must map a space list to itself (same spaces and dimensions in and out).');
end
% Two or more variables: the separable inverse (Woodbury)                   % MMP, 11/04/2026
if numel(P.vars)>1                                                          % MMP, 11/04/2026
    [Pinv,info] = inv_sep(P,opts);                                          % MMP, 11/04/2026
    return                                                                  % MMP, 11/04/2026
end                                                                         % MMP, 11/04/2026
isL2 = any(P.space_out,2);
l = find(isL2);     r = find(~isL2);
if numel(l)>1 || numel(r)>1
    error('copvar:inv:spaces','At most one R^k and one L2[s] space are supported; the grid is %d x %d.',M,Nc);
end
C = cell(M,M);
info = struct();
if isempty(l)
    % R^k -> R^k alone
    k = P.dim_out(r);
    C{r,r} = matblock(blockmat(P.C{r,r},k)\eye(k),k);
elseif isempty(r)
    [C{l,l},info] = inv(P.C{l,l},opts);
else
    k = P.dim_out(r);
    [Rh,info] = inv(P.C{l,l},opts);
    Q1 = P.C{r,l};  Q2 = P.C{l,r};
    T = blockmat(P.C{r,r},k);
    if ~isempty(Q1) && ~isempty(Q2)
        T = T - blockmat(Q1*Rh*Q2,k);
    end
    info.condT = cond(T);
    Th = matblock(T\eye(k),k);
    C{r,r} = Th;
    if ~isempty(Q1),    C{r,l} = -(Th*(Q1*Rh));     end
    if ~isempty(Q2),    C{l,r} = -((Rh*Q2)*Th);     end
    if ~isempty(Q1) && ~isempty(Q2)
        C{l,l} = Rh + (Rh*Q2)*(Th*(Q1*Rh));
    else
        C{l,l} = Rh;
    end
end
Pinv = copvar(C);
end


function Mat = blockmat(B,k)
% MAT = BLOCKMAT(B,K) the k x k matrix of an R^k -> R^k block, [] being zero.
% INPUT  B: 'sopvar' block on R^k, or []; k: its dimension.
% OUTPUT Mat: full k x k double.
if isempty(B),  Mat = zeros(k);     return,     end
Mat = full(B.params{1});
if ~isequal(size(Mat),[k k])
    error('copvar:inv:matblock','An R^%d -> R^%d block stores a %d x %d parameter.',k,k,size(Mat,1),size(Mat,2));
end
end

function B = matblock(Mat,k)
% B = MATBLOCK(MAT,K) the 'sopvar' block of a k x k matrix on R^k, as
% 'opvar2sopvar' builds it.
% INPUT  Mat: k x k double; k: its dimension.
% OUTPUT B: 'sopvar' block with no spatial variable.
vars = struct;  vars.in = {};   vars.out = {};
dom = struct;   dom.in = zeros(0,2);    dom.out = zeros(0,2);
B = sopvar({Mat},vars,cell(1,0),cell(1,0),dom,[k k]);
end


% BEGIN MMP, 11/04/2026: separable inverse in two or more variables (new)
function [Pinv,info] = inv_sep(P,opts)
% [PINV,INFO] = INV_SEP(P,OPTS) the inverse of a square separable 'copvar'
% P = M + L H R in two or more spatial variables, by the Woodbury identity
%   P^{-1} = M^{-1} + M^{-1} L Hh R M^{-1},  Hh = -H (I + K H)^{-1},
%   K = R M^{-1} L.
% INPUT  P: 'copvar', square, separable (see INV); opts: 'struct' with the
%        optional fields tol, deg0, degmax, deg, nq, septol (see INV).
% OUTPUT Pinv: 'copvar' over the spaces of P; info: 'struct' (see INV).

% Options and their defaults
o = struct('tol',1e-8,'deg0',4,'degmax',14,'nq',40,'septol',1e-10);
for f = fieldnames(o)'
    if isfield(opts,f{1}) && ~isempty(opts.(f{1})),     o.(f{1}) = opts.(f{1});     end
end
if isfield(opts,'deg') && ~isempty(opts.deg),   o.deg0 = opts.deg;  o.degmax = opts.deg;    end
% Each space: its variables (registry order), domain and dimension
nsp = size(P,1);    q = P.dim_out(:);
V = cell(nsp,1);    D = cell(nsp,1);
for i = 1:nsp
    V{i} = P.vars(P.space_out(i,:));    D{i} = P.dom(P.space_out(i,:),:);
end

% 1. Split P = M + F block by block, testing separability: the multiplier
% index (1 in every shared variable) only on the diagonal, mixed indices
% zero, and every integral index (2 or 3 in each shared variable) equal
Mb = cell(nsp,1);
F = struct('i',{},'j',{},'C',{},'ZL',{},'ZR',{},'vo',{},'vi',{});
for i = 1:nsp
    for j = 1:nsp
        B = P.C{i,j};
        if isempty(B),  continue,   end
        pr = B.params;
        scale = max(cellfun(@(c) full(max([abs(c(:));0])),pr(:)));
        if scale==0,    continue,   end
        k = numel(intersect(B.vars.out,B.vars.in));
        if k==0
            % No shared variable: the matrix of an R space, or a finite-rank
            % term between two spaces
            if i==j,    Mb{i} = blk(pr{1},B);   else,   F(end+1) = term(i,j,pr{1},B);  end %#ok<AGROW>
            continue
        end
        Cint = [];
        for n = 1:numel(pr)
            a = alpha(size(pr),n,k);    c = pr{n};
            if all(a==1)
                if i==j,    Mb{i} = blk(c,B);
                elseif any(abs(c(:))>o.septol*scale)
                    error('copvar:inv:separable','Block (%d,%d) multiplies between two spaces: P is not separable.',i,j);
                end
            elseif all(a>=2)
                if isempty(Cint),   Cint = c;
                elseif any(abs(c(:)-Cint(:))>o.septol*scale)
                    error('copvar:inv:separable','Block (%d,%d) has different lower and upper kernels: P is not separable.',i,j);
                end
            elseif any(abs(c(:))>o.septol*scale)
                error('copvar:inv:separable','Block (%d,%d) multiplies in one variable and integrates in another: P is not separable.',i,j);
            end
        end
        if ~isempty(Cint) && any(Cint(:)),  F(end+1) = term(i,j,Cint,B);  end %#ok<AGROW>
    end
end
for i = 1:nsp
    if isempty(Mb{i})
        error('copvar:inv:multiplier','Block (%d,%d) has no multiplier: the separable inverse needs an invertible M.',i,i);
    end
end
% H = blkdiag of the integral coefficients; row and column offsets of each
% term (rows: L factors, columns: R factors)
nF = numel(F);  nb = zeros(nF,1);   mb = zeros(nF,1);
for b = 1:nF,   [nb(b),mb(b)] = size(F(b).C);   end
roff = [0;cumsum(nb)];  coff = [0;cumsum(mb)];
H = zeros(sum(nb),sum(mb));
for b = 1:nF,   H(roff(b)+(1:nb(b)),coff(b)+(1:mb(b))) = F(b).C;    end

% 2. K = R M^{-1} L by quadrature, M^{-1} inverted at each node; repeated
% with 10 more nodes per variable to estimate the quadrature error
K = kmat(o.nq);     K2 = kmat(o.nq+10);
Kdiff = norm(K-K2,'fro')/max(norm(K2,'fro'),realmin);
% 3. Hh = -H (I + K H)^{-1}
IKH = eye(sum(mb)) + K*H;
Hh = -H/IKH;

% 4. Fit M^{-1} [I, L] and R M^{-1} on each space at degree d per
% variable, raised until the relative RMS residual is below tol
for d = o.deg0:2:max(o.deg0,o.degmax)
    [G,Lam,rr] = fits(d);
    if rr<=o.tol,   break,  end
end

% 5. The blocks of P^{-1}: the kernel Gam_i Hh_ij Lam_j in every integral
% index, and on the diagonal the multiplier M_i^{-1}
C = cell(nsp,nsp);
for i = 1:nsp
    Li = find([F.i]==i);    rows = cell2mat(arrayfun(@(b) roff(b)+(1:nb(b)),Li(:)','UniformOutput',false));
    for j = 1:nsp
        Rj = find([F.j]==j);    cols = cell2mat(arrayfun(@(b) coff(b)+(1:mb(b)),Rj(:)','UniformOutput',false));
        Cij = [];
        if ~isempty(rows) && ~isempty(cols)
            Cij = G{i}(:,q(i)+1:end)*Hh(rows,cols)*Lam{j};
        end
        % an off-diagonal block without a Woodbury term is zero
        if i~=j && (isempty(Cij) || ~any(Cij(:))),     continue,   end
        C{i,j} = pblock(i,j,Cij,d);
    end
end
Pinv = copvar(C);
info = struct('method','separable','d',d,'relrms',rr,'condIKH',cond(IKH),'Kdiff',Kdiff,'nF',nF, ...
              'note','Woodbury identity on M + L H R; M^{-1} [I, L] and R M^{-1} by tensor Chebyshev least squares');
if rr>o.tol
    info.note = sprintf('%s; fit tolerance %.1e NOT met at degmax %d (relrms %.2e)',info.note,o.tol,d,rr);
end

    function K = kmat(nq)
    % K = KMAT(NQ) the matrix R M^{-1} L (sum(mb) x sum(nb)): the block of
    % R factor b and L factor b', where b integrates over the space that b'
    % maps into, is the integral over that space of R_b M^{-1} L_b', at NQ
    % Gauss-Legendre nodes per variable.
    K = zeros(sum(mb),sum(nb));
    for jj = 1:nsp
        bs = find([F.j]==jj);   bp = find([F.i]==jj);
        if isempty(bs) || isempty(bp),  continue,   end
%       [X,w] = gl_grid(D{jj},nq);                                          % MMP, 11/09/2026 (was)
        % The tensor Gauss-Legendre rule of nq nodes per variable on the box % MMP, 11/09/2026
        % D{jj}: points X (one row each) and weights w; with no variable, one % MMP, 11/09/2026
        % point with weight 1 (the former locals GL_GRID and GAUSS1, the    % MMP, 11/09/2026
        % latter Golub-Welsch on [a,b] for each variable)                   % MMP, 11/09/2026
        if size(D{jj},1)==0                                                 % MMP, 11/09/2026
            X = zeros(1,0);     w = 1;                                      % MMP, 11/09/2026
        else                                                                % MMP, 11/09/2026
            g = cell(1,size(D{jj},1));    w = 1;                            % MMP, 11/09/2026
            for v = 1:size(D{jj},1)                                         % MMP, 11/09/2026
                k = (1:nq-1)';   bt = k./sqrt(4*k.^2-1);                    % MMP, 11/09/2026
                [Vv,Dv] = eig(diag(bt,1)+diag(bt,-1));                      % MMP, 11/09/2026
                [t,ix] = sort(diag(Dv));                                    % MMP, 11/09/2026
                g{v} = (D{jj}(v,1)+D{jj}(v,2))/2 + (D{jj}(v,2)-D{jj}(v,1))/2*t; % MMP, 11/09/2026
                wv = (D{jj}(v,2)-D{jj}(v,1))/2*2*Vv(1,ix)'.^2;              % MMP, 11/09/2026
                w = kron(w,wv);                                             % MMP, 11/09/2026
            end                                                             % MMP, 11/09/2026
            X = tensor_points(g);                                           % MMP, 11/09/2026
        end                                                                 % MMP, 11/09/2026
        Mi = minv_at(jj,X);
        ZRb = arrayfun(@(b) zeval(F(b).ZR,F(b).vi,V{jj},X),bs,'UniformOutput',false);
        ZLb = arrayfun(@(b) zeval(F(b).ZL,F(b).vo,V{jj},X),bp,'UniformOutput',false);
        for u = 1:numel(bs)
            ru = coff(bs(u))+(1:mb(bs(u)));
            for v = 1:numel(bp)
                cv = roff(bp(v))+(1:nb(bp(v)));
                Kuv = zeros(numel(ru),numel(cv));
                for p = 1:numel(w)
                    Rb = kron(eye(q(jj)),ZRb{u}(p,:).');
                    Lb = kron(eye(q(jj)),ZLb{v}(p,:));
                    Kuv = Kuv + w(p)*(Rb*Mi(:,:,p)*Lb);
                end
                K(ru,cv) = K(ru,cv) + Kuv;
            end
        end
    end
    end

    function Mi = minv_at(i,X)
    % MI = MINV_AT(I,X) the inverse of the multiplier M_I at the points X
    % (one row each, columns the variables V{I}), q(I) x q(I) x size(X,1).
    ZLm = zeval(Mb{i}.ZL,Mb{i}.vo,V{i},X);  ZRm = zeval(Mb{i}.ZR,Mb{i}.vi,V{i},X);
    Mi = zeros(q(i),q(i),size(X,1));
    for p = 1:size(X,1)
        Mp = kron(eye(q(i)),ZLm(p,:))*Mb{i}.C*kron(eye(q(i)),ZRm(p,:).');
        Mi(:,:,p) = Mp\eye(q(i));
    end
    end

    function [G,Lam,rr] = fits(d)
    % [G,LAM,RR] = FITS(D) for each space i the coefficients G{i}
    % ((q_i nW) x (q_i + n_i)) of M_i^{-1} [I, L_b ...] and LAM{i}
    % (m_i x (q_i nW)) of [R_b; ...] M_i^{-1} on the tensor monomials of
    % degree D per variable (nW = (D+1)^k_i, 1 on R; exact values on R),
    % and RR the largest relative RMS residual of the fits.
    G = cell(nsp,1);    Lam = cell(nsp,1);  rr = 0;
    for i = 1:nsp
        Li = find([F.i]==i);    Ri = find([F.j]==i);    kv = numel(V{i});
%       if kv==0,   X = zeros(1,0);     else,   X = cheb_grid(D{i},2*(d+1));    end % MMP, 11/09/2026 (was)
        % X: the tensor grid of 2(d+1) first-kind Chebyshev points per      % MMP, 11/09/2026
        % variable on the box D{i}, one row per point (the former local     % MMP, 11/09/2026
        % CHEB_GRID); no variable: one empty point                          % MMP, 11/09/2026
        if kv==0                                                            % MMP, 11/09/2026
            X = zeros(1,0);                                                 % MMP, 11/09/2026
        else                                                                % MMP, 11/09/2026
            nc = 2*(d+1);   g = cell(1,size(D{i},1));                       % MMP, 11/09/2026
            for v = 1:size(D{i},1)                                          % MMP, 11/09/2026
                g{v} = (D{i}(v,1)+D{i}(v,2))/2 + (D{i}(v,2)-D{i}(v,1))/2*cos(pi*((1:nc)'-0.5)/nc); % MMP, 11/09/2026
            end                                                             % MMP, 11/09/2026
            X = tensor_points(g);                                           % MMP, 11/09/2026
        end                                                                 % MMP, 11/09/2026
        Mi = minv_at(i,X);  np = size(X,1);
        ZLi = arrayfun(@(b) zeval(F(b).ZL,F(b).vo,V{i},X),Li,'UniformOutput',false);
        ZRi = arrayfun(@(b) zeval(F(b).ZR,F(b).vi,V{i},X),Ri,'UniformOutput',false);
        ncA = q(i) + sum(nb(Li));   nrB = sum(mb(Ri));
        A = zeros(np,q(i)*ncA);     Bv = zeros(np,nrB*q(i));
        for p = 1:np
            Lp = eye(q(i));
            for u = 1:numel(Li),    Lp = [Lp, kron(eye(q(i)),ZLi{u}(p,:))];     end %#ok<AGROW>
            Ap = Mi(:,:,p)*Lp;      A(p,:) = Ap(:).';
            if nrB>0
                Rp = zeros(0,q(i));
                for u = 1:numel(Ri),    Rp = [Rp; kron(eye(q(i)),ZRi{u}(p,:).')];   end %#ok<AGROW>
                Bp = Rp*Mi(:,:,p);  Bv(p,:) = Bp(:).';
            end
        end
        if kv==0
            % R: the values themselves, on the constant monomial
            cA = A;     cB = Bv;    nW = 1;
        else
            [cA,rA] = cheb_fit(X,A,D{i},d);     rr = max(rr,rA);
            cB = zeros((d+1)^kv,0);
            if nrB>0,   [cB,rB] = cheb_fit(X,Bv,D{i},d);    rr = max(rr,rB);    end
            nW = (d+1)^kv;
        end
        % entry (a,c) of the q x ncA matrix -> rows (a-1)nW+(1:nW) of column c
        G{i} = zeros(q(i)*nW,ncA);
        for c = 1:ncA
            for a = 1:q(i),     G{i}((a-1)*nW+(1:nW),c) = cA(:,a+(c-1)*q(i));    end
        end
        % entry (c,a) of the nrB x q matrix -> columns (a-1)nW+(1:nW) of row c
        Lam{i} = zeros(nrB,q(i)*nW);
        for a = 1:q(i)
            for c = 1:nrB,      Lam{i}(c,(a-1)*nW+(1:nW)) = cB(:,c+(a-1)*nrB).';    end
        end
    end
    end

    function Bk = pblock(i,j,Cij,d)
    % BK = PBLOCK(I,J,CIJ,D) the 'sopvar' block (I,J) of P^{-1}: CIJ in every
    % integral index on the tensor monomials of degree D of the variables of
    % spaces I and J, and on the diagonal the multiplier M_I^{-1}, its
    % coefficients in the columns of the constant ZR monomial (canonical
    % multiplier form).
    ki = numel(V{i});   kj = numel(V{j});
    nWi = (d+1)^ki;     nWj = (d+1)^kj;
    if isempty(Cij),    Cij = zeros(q(i)*nWi,q(j)*nWj);     end
    C0 = zeros(q(i)*nWi,q(j)*nWj);
    if i==j
        Gm = G{i}(:,1:q(i));
        for a = 1:q(i),     C0(:,(a-1)*nWj+1) = Gm(:,a);    end
    end
    k = numel(intersect(V{i},V{j}));
    if k==0
        pr = {sparse(C0+Cij)};
    else
        pr = cell([3*ones(1,k),1]);
        for e = 1:numel(pr)
            a = alpha(size(pr),e,k);
            if all(a==1),       pr{e} = sparse(C0);
            elseif all(a>=2),   pr{e} = sparse(Cij);
            else,               pr{e} = sparse(size(C0,1),size(C0,2));
            end
        end
    end
    vars = struct;  vars.out = V{i};    vars.in = V{j};
    dom = struct;   dom.out = D{i};     dom.in = D{j};
    Bk = sopvar(pr,vars,repmat({(0:d)'},1,ki),repmat({(0:d)'},1,kj),dom,[q(i) q(j)]);
    end
end


function s = blk(C,B)
% S = BLK(C,B) the multiplier parameter C of the 'sopvar' block B, with
% the block's monomial exponents ZL, ZR and its variables vo (out), vi (in).
s = struct('C',full(C),'ZL',{B.ZL},'ZR',{B.ZR},'vo',{cellstr(B.vars.out)},'vi',{cellstr(B.vars.in)});
end

function t = term(i,j,C,B)
% T = TERM(I,J,C,B) the integral term of block (I,J): the coefficient C of
% the kernel (I x Z_L)' C (I x Z_R) of the 'sopvar' block B, with its
% exponents and variables, as a 'struct' with fields i, j, C, ZL, ZR, vo, vi.
t = struct('i',i,'j',j,'C',full(C),'ZL',{B.ZL},'ZR',{B.ZR},'vo',{cellstr(B.vars.out)},'vi',{cellstr(B.vars.in)});
end

function a = alpha(sz,n,k)
% A = ALPHA(SZ,N,K) the multi-index (1 multiplier, 2 lower, 3 upper in each
% shared variable) of the N-th entry of a params cell of size SZ with K
% shared variables.
if k==1,    a = n;  return,     end
c = cell(1,numel(sz));  [c{:}] = ind2sub(sz,n);     a = [c{1:k}];
end

function Z = zeval(Zc,bv,Vs,X)
% Z = ZEVAL(ZC,BV,VS,X) the values of the monomial vector kron_v
% BV{v}.^ZC{v} at the points X (one row each, columns the variables VS),
% one row per point; a column of ones for no variable.
Z = ones(size(X,1),1);
for v = 1:numel(bv)
    Zv = X(:,strcmp(Vs,bv{v})).^(Zc{v}(:).');
    n1 = size(Z,2);     n2 = size(Zv,2);
    Z = Z(:,kron(1:n1,ones(1,n2))).*Zv(:,repmat(1:n2,1,n1));
end
end

function [cM,rr] = cheb_fit(X,Fv,dom,d)
% [CM,RR] = CHEB_FIT(X,FV,DOM,D) least squares fit of the columns of FV
% (values at the points X, one column of X per variable, box DOM) by the
% tensor Chebyshev polynomials of degree D per variable, returned as the
% coefficients CM of the tensor monomials kron_v x_v.^(0:D)', and RR, the
% relative RMS residual at X.
Phi = ones(size(X,1),1);    T = 1;
for v = 1:size(X,2)
%   [Pv,Tv] = cheb1(X(:,v),dom(v,:),d);                                     % MMP, 11/09/2026 (was)
    % Pv: the values at X(:,v) of the Chebyshev polynomials T_0..T_d of     % MMP, 11/09/2026
    % t = (2x-a-b)/(b-a) on dom(v,:) = [a b]; Tv: column k+1 the monomial   % MMP, 11/09/2026
    % coefficients of T_k in x, ascending degree (the former local CHEB1)   % MMP, 11/09/2026
    xv = X(:,v);                                                            % MMP, 11/09/2026
    al = 2/(dom(v,2)-dom(v,1));     be = -(dom(v,1)+dom(v,2))/(dom(v,2)-dom(v,1)); % MMP, 11/09/2026
    t = al*xv(:)+be;                                                        % MMP, 11/09/2026
    Pv = ones(numel(xv),d+1);   Tv = zeros(d+1,d+1);    Tv(1,1) = 1;        % MMP, 11/09/2026
    if d>=1,    Pv(:,2) = t;    Tv(1:2,2) = [be;al];    end                 % MMP, 11/09/2026
    for k = 2:d                                                             % MMP, 11/09/2026
        Pv(:,k+1) = 2*t.*Pv(:,k) - Pv(:,k-1);                               % MMP, 11/09/2026
        Tv(:,k+1) = 2*(al*[0;Tv(1:end-1,k)] + be*Tv(:,k)) - Tv(:,k-1);      % MMP, 11/09/2026
    end                                                                     % MMP, 11/09/2026
    n1 = size(Phi,2);   n2 = d+1;
    Phi = Phi(:,kron(1:n1,ones(1,n2))).*Pv(:,repmat(1:n2,1,n1));
    T = kron(T,Tv);
end
cP = Phi\Fv;
cM = T*cP;
rr = norm(Phi*cP-Fv,'fro')/max(norm(Fv,'fro'),realmin);
end

function X = tensor_points(g)
% X = TENSOR_POINTS(G) all combinations of the node vectors in the cell G,
% one row per point, the last variable running fastest (the order of kron).
X = zeros(1,0);
for v = 1:numel(g)
    n1 = size(X,1);     n2 = numel(g{v});
    X = [X(kron(1:n1,ones(1,n2)),:), reshape(g{v}(repmat(1:n2,1,n1)),[],1)];
end
end

% END MMP, 11/04/2026
