function W = pie_witness_sop(PIE,mode,opts,Y,wblk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W = PIE_WITNESS_SOP(PIE,MODE,OPTS,Y,WBLK) computes, for the 1-D PIE
%
%       T x_t + Tw w_t = A x + B1 w,      z = C1 x + D11 w,
%
% the quantity MODE names on the Galerkin matrices T, Tw, A, B, C, D of its
% operators in the orthonormal Legendre basis of degree N = OPTS.N_cheb
% (PIE_DISC_SOP), in which the Euclidean norm of a coefficient vector is
% the norm of the function it represents:
%
%   'gain'  max norm(G(om)), G(om) = C*((1i*om*T - A)\(B - 1i*om*Tw)) + D,
%           over om = 0
%           and logspace(-3,log10(100*rho),OPTS.nfreq), rho the largest of
%           1 and the moduli of the finite eigenvalues of (A,T), refined by
%           fminbnd between the two samples next to the largest; the worst
%           input w0 is the right singular vector of G at the maximizing
%           frequency;
%   'rate'  the largest real part of the finite eigenvalues of (A,T);
%   'h2'    sqrt(trace(C*X*C')), X = lyap(T\A,(T\B)*(T\B)'), D not used;
%           computed only when Tw is zero.
%
% In MODE 'gain' with Y, it also returns the alignment of the diagonal
% block Yb = Y.C{WBLK,WBLK} with w0,
%
%   abs(sum(sum(K.*S)))/(norm(K,'fro')*norm(S,'fro')),   S = real(w0*w0'),
%
% with K = (Yb + Yb')/2 and w0 the coefficient vector when w is finite;
% when w is distributed, K the kernel of Yb and S at the N+1 Gauss-Legendre
% nodes of the domain, the multiplier of Yb on the diagonal divided by h_i,
% and both weighted entrywise by sqrt(h_i*h_j), h_i the quadrature weight
% of node i.
%
% INPUT
% - PIE:    'pie_struct' object of a 1-D PIE with T, Tw, A, B1, C1, D11
%           'opvar' or 'copvar' objects (the other operators are not used);
%           w has nwf finite and nwd distributed components,
%           PIE.B1.dim(:,2) after conversion to 'opvar' objects;
% - mode:   char, 'gain', 'rate' or 'h2';
% - opts:   (optional) 'struct' with the optional fields N_cheb (double,
%           the Legendre degree N, default 24) and nfreq (double, default
%           240);
% - Y:      (optional) 'copvar' object whose diagonal block (WBLK,WBLK) is
%           a 'sopvar' object on the space of w: on R^nwf its matrix, on
%           L2^nwd[s] its multiplier, lower kernel and, if stored, upper
%           kernel (else the adjoint of the lower kernel); default []: no
%           alignment;
% - wblk:   (optional) double, the index of that block (default 1).
%
% OUTPUT
% - W:      'struct' with fields
%             ok          logical, true when the quantity of MODE was
%                         computed;
%             note        char, the reason when ok is false (a 2-D PIE, a
%                         failed discretization or Lyapunov solve, Tw
%                         nonzero in 'h2') and the message of a failed
%                         alignment;
%             N           double, the Legendre degree (NaN without a
%                         discretization);
%             gain, omega double, the largest sampled norm of G and its
%                         frequency ('gain');
%             wcoef       complex double column, w0: the nwf finite
%                         entries, then the N+1 Legendre coefficients of
%                         each of the nwd distributed components ('gain');
%             sgrid       double column, the N+1 Gauss-Legendre nodes of
%                         the domain ('gain');
%             wgrid       (N+1) x nwd complex double, the distributed
%                         components of w0 at sgrid ('gain'; [] when
%                         nwd = 0);
%             maxre, rate double, the largest real part and -maxre
%                         ('rate');
%             h2          double ('h2');
%             dual_alignment  double in [0,1] ('gain' with Y; NaN when w
%                         has finite and distributed components or the
%                         block is empty);
%             t           double, the computation time in seconds;
%           a field that is not computed is NaN or [].
%
% NOTES
% For a 1-D PIE, a MODE other than 'gain', 'rate' and 'h2' is an error.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_witness_sop
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
% Initial coding MMP, 10/08/2026
% MMP, 10/28/2026: Header to the header protocol (what the function does,
%                inputs, outputs), the licence block added, step comments
%                in the body and headers of the local functions; code
%                unchanged. The deleted header text was the interpretation
%                of the outputs (a lower bound, the meaning of alignment
%                values) and the reason the discretization is exact.
% MMP, 10/28/2026: 'gain' includes the disturbance at the boundary: from
%                T x_t + Tw w_t = A x + B1 w, the input operator of G(om) is
%                B - i*om*Tw. Before, Tw was not used, so for Tw ~= 0 the
%                witness was the gain of the PIE without the Tw term.
%                'h2' returns ok = false when Tw ~= 0. For Tw = 0 both are
%                computed as before.
% MMP, 10/28/2026: 'gain' with an empty channel (no disturbance or no
%                regulated output) returns the gain 0 at omega 0. Before,
%                the worst input was read off the SVD of an nz x 0 response
%                and the index failed, so the H-infinity gain executive
%                errored on a PIE without disturbance (maintainer: empty
%                channels should not raise an error).
% DJ, 10/29/2026: Convert the operators of the PIE to 'opvar' objects
%                using pie2opvar, which PIE_DISC_SOP and PIE.B1.dim require,
%                instead of initializing the PIE, whose operators are now
%                'copvar' objects. This is a temporary fix, we should 
%                really work with the new PIE structure directly;

% Output with every quantity not computed
W = struct('ok',false,'note','','N',NaN,'gain',NaN,'omega',NaN,'wcoef',[],'wgrid',[],'sgrid',[], ...
           'maxre',NaN,'rate',NaN,'h2',NaN,'dual_alignment',NaN,'t',0);
% Defaults of the optional inputs; the degree and frequency count of OPTS
if nargin<3 || isempty(opts),   opts = struct();    end
if nargin<4,    Y = [];     end
if nargin<5,    wblk = 1;   end
N = getf(opts,'N_cheb',24);     nfreq = getf(opts,'nfreq',240);
% Start the timer; convert the PIE to opvar objects (an error leaves it as
% given)
t0 = tic;
try
    PIE = pie2opvar(PIE);                                                   % DJ, 10/29/2026
catch
end
% 1-D only
if PIE.dim~=1
    W.note = '2-D PIE: no discretization in this version';     return
end
% Galerkin matrices of T, A, B1, C1, D11; D zero when D11 is absent or
% does not match C*B
try
    [T,di] = pie_disc_sop(PIE.T,N);
    A = pie_disc_sop(PIE.A,N);
    B = pie_disc_sop(PIE.B1,N);
    C = pie_disc_sop(PIE.C1,N);
    D = pie_disc_sop(PIE.D11,N);
    % The boundary-disturbance operator Tw, when nonzero                    % MMP, 10/28/2026
    Tw0 = (PIE.Tw==0);                                                      % MMP, 10/28/2026
    Twd = [];   if ~Tw0,    Twd = pie_disc_sop(PIE.Tw,N);   end             % MMP, 10/28/2026
catch ME
    W.note = ['pie_disc_sop: ' ME.message];     return
end
if isempty(D) || ~isequal(size(D),[size(C,1) size(B,2)]),  D = zeros(size(C,1),size(B,2));    end
W.N = N;
% Domain and number of basis functions; finite and distributed components
% of w
a = di.a;   b = di.b;   q = di.q;
nwf = PIE.B1.dim(1,2);  nwd = PIE.B1.dim(2,2);
nstate = size(A,1);
switch mode
    case 'rate'
        % Largest real part of the finite eigenvalues of (A,T)
        lam = eig(A,T);
        lam = lam(isfinite(lam));
        W.maxre = max(real(lam));   W.rate = -W.maxre;    W.ok = true;
    case 'h2'
        % No H2 witness for a disturbance at the boundary (Tw ~= 0)         % MMP, 10/28/2026
        if ~Tw0                                                             % MMP, 10/28/2026
            W.note = 'h2: not computed for Tw ~= 0';                        % MMP, 10/28/2026
            W.t = toc(t0);  return                                          % MMP, 10/28/2026
        end                                                                 % MMP, 10/28/2026
        % Lyapunov equation of (T\A, T\B), then the trace with C
        try
            At = T\A;   Bt = T\B;
            X = lyap(At,Bt*Bt');
            W.h2 = sqrt(abs(trace(C*X*C')));   W.ok = true;
        catch ME
            W.note = ['lyap: ' ME.message];
        end
    case 'gain'
        % No disturbance or no regulated output: the response is empty and  % MMP, 10/28/2026
        % the gain 0                                                        % MMP, 10/28/2026
        if size(B,2)==0 || size(C,1)==0                                     % MMP, 10/28/2026
            W.omega = 0;    W.gain = 0;     W.ok = true;                    % MMP, 10/28/2026
            W.t = toc(t0);  return                                          % MMP, 10/28/2026
        end                                                                 % MMP, 10/28/2026
        % Norm of the frequency response at om = 0 and nfreq frequencies up
        % to 100*rho, rho the largest of 1 and the finite eigenvalue moduli
        G = @(om) C*((1i*om*T - A)\B) + D;
        % With Tw ~= 0 the input operator is B - i*om*Tw                    % MMP, 10/28/2026
        if ~Tw0,    G = @(om) C*((1i*om*T - A)\(B - 1i*om*Twd)) + D;    end % MMP, 10/28/2026
        lam = eig(A,T);     lam = lam(isfinite(lam));
        rho = max([abs(lam); 1]);
        og = [0, logspace(-3,log10(100*rho),nfreq)];
        sg = zeros(size(og));
        for k = 1:numel(og),    sg(k) = norm(G(og(k)));    end
        % Refine the largest sample between its two neighbours
        [~,k] = max(sg);
        lo = og(max(k-1,1));    hi = og(min(k+1,numel(og)));
        W.omega = og(k);    W.gain = sg(k);
        if hi>lo
            f = @(x) -norm(G(x));
            [om,sv] = fminbnd(f,lo,hi,optimset('TolX',1e-7*max(hi,1),'Display','off'));
            if -sv>=sg(k),  W.omega = om;   W.gain = -sv;   end
        end
        % Worst input: the right singular vector at the maximizing frequency
        [~,~,V] = svd(G(W.omega));
        wc = V(:,1);                                % worst input, orthonormal coefficients
        W.wcoef = wc;
        % Its distributed components at the Gauss-Legendre nodes
        [sq,wq] = gauss_legendre(q,a,b);
        W.sgrid = sq;
        if nwd>0
%           Pq = (sq(:).^(0:N))*di.C';              % p_j(s_q): q x (N+1)  % MMP, 10/08/2026 (was)
            Pq = di.evalbasis(sq(:));               % p_j(s_q) by recurrence % MMP, 10/08/2026
            wd = reshape(wc(nwf+1:end),q,nwd);
            W.wgrid = Pq*wd;                        % values, one column per component
        end
        W.ok = true;
        % Alignment of block WBLK of Y with the worst input
        if ~isempty(Y)
            try
                W.dual_alignment = dual_alignment(Y,wblk,W,nwf,nwd,sq,wq);
            catch ME
                W.note = [W.note ' | dual alignment: ' ME.message];
            end
        end
    otherwise
        error('pie_witness_sop:mode','MODE should be ''gain'', ''rate'' or ''h2''.');
end
W.t = toc(t0);
end


function v = getf(s,f,d)
% V = GETF(S,F,D) returns the field F of the struct S, or D when S has no
% nonempty field F.
% INPUT  s: 'struct'; f: char, the field name; d: the default value.
% OUTPUT v: the field's value or D.
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = s.(f);  else,   v = d;  end
end


function [x,w] = gauss_legendre(n,a,b)
% [X,W] = GAUSS_LEGENDRE(N,A,B) returns the N-point Gauss-Legendre nodes and
% weights on [A,B], by the Golub-Welsch method.
% INPUT  n: double, the number of nodes; a, b: double, the interval.
% OUTPUT x: n x 1 double, the nodes in increasing order;
%        w: 1 x n double, the weights.
% Nodes: the eigenvalues of the Jacobi matrix of the Legendre recurrence;
% weights: twice the squared first entries of its eigenvectors; then the
% map from [-1,1] to [a,b]
beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
J = diag(beta,1)+diag(beta,-1);
[V,Dg] = eig(J);
[x,i] = sort(diag(Dg));     w = 2*V(1,i).^2;
x = (b-a)/2*x + (a+b)/2;    w = (b-a)/2*w;
end


function al = dual_alignment(Y,wblk,W,nwf,nwd,sq,wq)
% AL = DUAL_ALIGNMENT(Y,WBLK,W,NWF,NWD,SQ,WQ) returns the normalized
% Frobenius inner product of the symmetrized block (WBLK,WBLK) of Y with
% the covariance real(w0*w0') of the worst input w0: on R^NWF, or at the
% nodes SQ with each entry weighted by the quadrature weights WQ.
% INPUT  Y:        'copvar' object; wblk: double, the block index;
%        W:        'struct' of PIE_WITNESS_SOP, fields wcoef and wgrid;
%        nwf, nwd: double, the finite and distributed components of w;
%        sq, wq:   double, the Gauss-Legendre nodes and weights.
% OUTPUT al:       double in [0,1]; NaN for an empty block or for w with
%                  finite and distributed components.
% The block; no comparison when it is empty or w is mixed
Yb = Y.C{wblk,wblk};
if isempty(Yb),     al = NaN;   return,     end
if nwf>0 && nwd>0,  al = NaN;   return,     end     % mixed space: not compared
if nwd==0                                           % finite disturbance: the multiplier matrix
    K = full(block_matrix(Yb));
    K = (K+K')/2;
    wf = W.wcoef(1:nwf);
    S = real(wf*wf');
    al = abs(sum(sum(K.*S)))/max(norm(K,'fro')*norm(S,'fro'),eps);
    return
end
% Distributed w: kernel and covariance at the nodes, entry (i,j) weighted
% by sqrt(wq_i*wq_j)
%Kq = kernel_on_grid(Yb,sq);                        % nwd q x nwd q, symmetrised  % MMP, 10/08/2026 (was)
Kq = kernel_on_grid(Yb,sq,wq);                      % nwd q x nwd q, symmetrised  % MMP, 10/08/2026
S = real(W.wgrid(:)*W.wgrid(:)');
sw = sqrt(repmat(wq(:),nwd,1));
Kw = (sw*sw').*Kq;  Sw = (sw*sw').*S;
al = abs(sum(sum(Kw.*Sw)))/max(norm(Kw,'fro')*norm(Sw,'fro'),eps);
end


function Kmat = block_matrix(B)
% KMAT = BLOCK_MATRIX(B) returns the matrix of a 'sopvar' block between
% finite-dimensional spaces: its first parameter, or zeros(B.dims) when
% that is empty.
% INPUT  B:    'sopvar' object without spatial variables.
% OUTPUT Kmat: B.dims(1) x B.dims(2) double.
p = B.params;
if iscell(p),   Kmat = p{1};    else,   Kmat = p;   end
if isempty(Kmat),   Kmat = zeros(B.dims);   end
end


function Kq = kernel_on_grid(B,sq,wq)                                       % MMP, 10/08/2026
% KQ = KERNEL_ON_GRID(B,SQ,WQ) returns the kernel of the 1-D fixed sopvar
% block B on L2^n[s] at the grid points SQ, symmetrized:
% lower cell for t < s, the stored upper cell or the lower cell's adjoint for
% t > s; the multiplier cell, a delta, goes on the diagonal divided by the  % MMP, 10/08/2026
% quadrature weight, so that the weighted Frobenius product below         % MMP, 10/08/2026
% integrates it as int y0(s) |w(s)|^2 ds. Component-outer order.          % MMP, 10/08/2026
% INPUT  B:   'sopvar' object, L2^n[s] -> L2^n[s], with parameters
%             {C0, C1} or {C0, C1, C2} (multiplier, t < s, t > s) and
%             monomial bases ZL{1}, ZR{1} in s;
%        sq:  double, the q grid points;
%        wq:  double, their q quadrature weights.
% OUTPUT Kq:  n*q x n*q double, symmetric, component-outer.
% Sizes, and the monomials of the bases at the grid points
p = B.params;   n = B.dims(1);
ZL = B.ZL{1}(:);    ZR = B.ZR{1}(:);
NL = numel(ZL);     NR = numel(ZR);
q = numel(sq);
zl = sq(:).^(ZL');  zr = sq(:).^(ZR');
% Multiplier, lower and upper cells; an empty or zero cell as a sparse zero
C0 = p{1};  if isempty(C0) || ~any(C0(:)),  C0 = sparse(n*NL,n*NR);     end % MMP, 10/08/2026
C1 = p{2};  C2 = [];    if numel(p)>=3,  C2 = p{3};    end
if isempty(C1) || ~any(C1(:)),  C1 = sparse(n*NL,n*NR);     end
% The n x n kernel block of every pair of grid points (i,j), grid-outer
Kg = zeros(n*q);                                    % grid-outer
for i = 1:q
    Li = kron(eye(n),zl(i,:));
    for j = 1:q
        Rj = kron(eye(n),zr(j,:)');
        % t < s: lower cell; t > s: upper cell, or the adjoint of the lower
        % cell; t = s: the lower cell symmetrized, plus the multiplier over
        % the quadrature weight
        if j<i
            K = Li*C1*Rj;
        elseif j>i
            if ~isempty(C2) && any(C2(:)),  K = Li*C2*Rj;
            else,   K = (kron(eye(n),zl(j,:))*C1*kron(eye(n),zr(i,:)'))';
            end
        else
            K = (Li*C1*Rj + (Li*C1*Rj)')/2;
            K = K + (Li*C0*Rj + (Li*C0*Rj)')/(2*wq(i));                     % MMP, 10/08/2026
        end
        Kg((i-1)*n+(1:n),(j-1)*n+(1:n)) = full(K);
    end
end
% Reorder component-outer and symmetrize
P = reshape(1:n*q,n,q)';    P = P(:);
Kq = Kg(P,P);
Kq = (Kq+Kq')/2;
end
