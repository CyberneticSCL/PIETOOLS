function W = pie_witness_sop(PIE,mode,opts,Y,wblk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% W = PIE_WITNESS_SOP(PIE,MODE,OPTS[,Y,WBLK]) the numerical counterpart of a
% certificate, on the Galerkin discretization of the 1-D PIE in the
% orthonormal Legendre basis of degree OPTS.N_cheb (PIE_DISC_SOP; the
% kernels are polynomial, so the matrices are exact and the Euclidean norm
% of a coefficient vector is the L2 norm):
%
%   MODE 'gain'  the frequency response G(iw) = C (iw T - A)^{-1} B + D, its
%                largest singular value over a frequency grid of OPTS.nfreq
%                points refined around the maximum: W.gain, W.omega, the
%                worst input W.wcoef (coefficients, complex) and its values
%                W.wgrid on the Gauss grid W.sgrid. A numerical LOWER bound
%                on the true gain, to set beside the certified upper bound;
%   MODE 'rate'  the generalized eigenvalues of (A,T): W.maxre the largest
%                real part, W.rate = -W.maxre the numerical decay rate;
%   MODE 'h2'    the H2 norm from the Lyapunov equation of (T^{-1}A, T^{-1}B,
%                C): W.h2.
%
% With Y (the dual kernel of the negativity constraint, LPIGETDUAL_SOP) and
% WBLK (the index of the disturbance space in that container) in MODE
% 'gain', W.dual_alignment is the normalised Frobenius inner product, on the
% grid, of the symmetrised dual kernel of block (WBLK,WBLK) with the
% covariance Re(w w^H) of the worst input: 1 when the solver's dual points
% at the worst input, 0 when orthogonal (distributed-only or finite-only
% disturbance spaces; NaN otherwise).
%
% W.ok is false, with W.note, for a 2-D PIE (no discretization in this
% version) or when the discretization fails.
%
% Initial coding MMP, 10/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W = struct('ok',false,'note','','N',NaN,'gain',NaN,'omega',NaN,'wcoef',[],'wgrid',[],'sgrid',[], ...
           'maxre',NaN,'rate',NaN,'h2',NaN,'dual_alignment',NaN,'t',0);
if nargin<3 || isempty(opts),   opts = struct();    end
if nargin<4,    Y = [];     end
if nargin<5,    wblk = 1;   end
N = getf(opts,'N_cheb',24);     nfreq = getf(opts,'nfreq',240);
t0 = tic;
try
    PIE = initialize(PIE);
catch
end
if PIE.dim~=1
    W.note = '2-D PIE: no discretization in this version';     return
end
try
    [T,di] = pie_disc_sop(PIE.T,N);
    A = pie_disc_sop(PIE.A,N);
    B = pie_disc_sop(PIE.B1,N);
    C = pie_disc_sop(PIE.C1,N);
    D = pie_disc_sop(PIE.D11,N);
catch ME
    W.note = ['pie_disc_sop: ' ME.message];     return
end
if isempty(D) || ~isequal(size(D),[size(C,1) size(B,2)]),  D = zeros(size(C,1),size(B,2));    end
W.N = N;
a = di.a;   b = di.b;   q = di.q;
nwf = PIE.B1.dim(1,2);  nwd = PIE.B1.dim(2,2);
nstate = size(A,1);
switch mode
    case 'rate'
        lam = eig(A,T);
        lam = lam(isfinite(lam));
        W.maxre = max(real(lam));   W.rate = -W.maxre;    W.ok = true;
    case 'h2'
        try
            At = T\A;   Bt = T\B;
            X = lyap(At,Bt*Bt');
            W.h2 = sqrt(abs(trace(C*X*C')));   W.ok = true;
        catch ME
            W.note = ['lyap: ' ME.message];
        end
    case 'gain'
        G = @(om) C*((1i*om*T - A)\B) + D;
        lam = eig(A,T);     lam = lam(isfinite(lam));
        rho = max([abs(lam); 1]);
        og = [0, logspace(-3,log10(100*rho),nfreq)];
        sg = zeros(size(og));
        for k = 1:numel(og),    sg(k) = norm(G(og(k)));    end
        [~,k] = max(sg);
        lo = og(max(k-1,1));    hi = og(min(k+1,numel(og)));
        W.omega = og(k);    W.gain = sg(k);
        if hi>lo
            f = @(x) -norm(G(x));
            [om,sv] = fminbnd(f,lo,hi,optimset('TolX',1e-7*max(hi,1),'Display','off'));
            if -sv>=sg(k),  W.omega = om;   W.gain = -sv;   end
        end
        [~,~,V] = svd(G(W.omega));
        wc = V(:,1);                                % worst input, orthonormal coefficients
        W.wcoef = wc;
        [sq,wq] = gauss_legendre(q,a,b);
        W.sgrid = sq;
        if nwd>0
%           Pq = (sq(:).^(0:N))*di.C';              % p_j(s_q): q x (N+1)  % MMP, 10/08/2026 (was)
            Pq = di.evalbasis(sq(:));               % p_j(s_q) by recurrence % MMP, 10/08/2026
            wd = reshape(wc(nwf+1:end),q,nwd);
            W.wgrid = Pq*wd;                        % values, one column per component
        end
        W.ok = true;
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
if isstruct(s) && isfield(s,f) && ~isempty(s.(f)),  v = s.(f);  else,   v = d;  end
end


function [x,w] = gauss_legendre(n,a,b)
beta = 0.5./sqrt(1-(2*(1:n-1)).^(-2));
J = diag(beta,1)+diag(beta,-1);
[V,Dg] = eig(J);
[x,i] = sort(diag(Dg));     w = 2*V(1,i).^2;
x = (b-a)/2*x + (a+b)/2;    w = (b-a)/2*w;
end


function al = dual_alignment(Y,wblk,W,nwf,nwd,sq,wq)
% Frobenius alignment between the symmetrised dual kernel of block (wblk,wblk)
% and the covariance of the worst input, both on the grid (or the finite space).
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
%Kq = kernel_on_grid(Yb,sq);                        % nwd q x nwd q, symmetrised  % MMP, 10/08/2026 (was)
Kq = kernel_on_grid(Yb,sq,wq);                      % nwd q x nwd q, symmetrised  % MMP, 10/08/2026
S = real(W.wgrid(:)*W.wgrid(:)');
sw = sqrt(repmat(wq(:),nwd,1));
Kw = (sw*sw').*Kq;  Sw = (sw*sw').*S;
al = abs(sum(sum(Kw.*Sw)))/max(norm(Kw,'fro')*norm(Sw,'fro'),eps);
end


function Kmat = block_matrix(B)
p = B.params;
if iscell(p),   Kmat = p{1};    else,   Kmat = p;   end
if isempty(Kmat),   Kmat = zeros(B.dims);   end
end


function Kq = kernel_on_grid(B,sq,wq)                                       % MMP, 10/08/2026
% The integral kernel of a 1-D fixed sopvar block on L2^n[s] on the grid:
% lower cell for t < s, the stored upper cell or the lower cell's adjoint for
% t > s; the multiplier cell, a delta, goes on the diagonal divided by the  % MMP, 10/08/2026
% quadrature weight, so that the weighted Frobenius product below         % MMP, 10/08/2026
% integrates it as int y0(s) |w(s)|^2 ds. Component-outer order.          % MMP, 10/08/2026
p = B.params;   n = B.dims(1);
ZL = B.ZL{1}(:);    ZR = B.ZR{1}(:);
NL = numel(ZL);     NR = numel(ZR);
q = numel(sq);
zl = sq(:).^(ZL');  zr = sq(:).^(ZR');
C0 = p{1};  if isempty(C0) || ~any(C0(:)),  C0 = sparse(n*NL,n*NR);     end % MMP, 10/08/2026
C1 = p{2};  C2 = [];    if numel(p)>=3,  C2 = p{3};    end
if isempty(C1) || ~any(C1(:)),  C1 = sparse(n*NL,n*NR);     end
Kg = zeros(n*q);                                    % grid-outer
for i = 1:q
    Li = kron(eye(n),zl(i,:));
    for j = 1:q
        Rj = kron(eye(n),zr(j,:)');
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
P = reshape(1:n*q,n,q)';    P = P(:);
Kq = Kg(P,P);
Kq = (Kq+Kq')/2;
end
