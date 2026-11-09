function Rd = bl_reduce(D,ref,opts)                                          % CC, 09/27/2026
% BL_REDUCE  Facial reduction of a SeDuMi-form SDP {At'x = b, x in K} and the
% map back.  The face comes from CERTIFICATES, not from a primal point.
%
%   Rd = bl_reduce(D,ref)          % D.At (nvar x m), D.b, D.K;  ref.y = a dual vector
%   Rd = bl_reduce(D,ref,opts)     % opts.sig (0.05), opts.drop (1e-9), opts.outmat
%
% FACE, two stages per PSD block:
%  1. exact coordinate zeros (bl_face): rows with b_i = 0 whose terms are
%     same-sign diagonals force those diagonals, and their rows/columns, to 0.
%     The basis is an exact selection, so those rows reduce to exactly zero.
%  2. the dual slack Z = -A'y (sign chosen so Z is ~PSD): on the remaining
%     coordinates, an eigen-direction v of Z_k is exposed and removed when its
%     Z-eigenvalue exceeds sig x the largest |eig(Z)| over ALL blocks AND (when
%     ref.x is given) the primal point is small along it, v'X_k v <= taux x
%     lambda_max(X_k) (complementarity).  b'y ~ 0 makes <X,Z> ~ 0 for every
%     feasible X, which forces X to vanish only where Z is large in absolute
%     terms.  Measured: a per-block relative threshold exposed directions of
%     blocks whose Z is ~1e-6..1e-9 overall (numerically zero), removed real
%     directions, and Mosek then certified the reduced program INFEASIBLE
%     (hinf_rd1, hinf_rd1_hv, hinfco_rd1).  Exposed directions have Z 0.3-1.0 of
%     the block max where X is ~1e-8; genuine ones have Z <= 1e-3.
% V_k = [selection] * [kept eigenvectors of Z].  Blocks with nothing to remove
% are left unchanged (no work).
%
% WHY NOT THE PRIMAL POINT (measured).  Taking V from X's large eigenvectors got
% the face dimension right but mixed the exposed coordinate directions into V at
% ~1e-2, so a face row "X_rr = 0" became <v v', Y> = 0 with v ~ 1e-2: again no
% interior, and Mosek returned UNKNOWN on the reduced program although Mosek's
% own projected point was strictly feasible for it (lambda_min(Y0) ~ 1e-5).
% CC, 09/27/2026 (b) CORRECTION: a projected Y0 = V'XV is PD but NOT feasible
% (row residual ~1e-5, the square root of the dropped ~1e-8 diagonals through
% the off-diagonals), and the min-norm repair turns it indefinite; "strictly
% feasible" above is unproven.  Test an interior with the margin program
% max t s.t. Y_k >= t s_k I instead.
%
% STATUS (CC, 09/27/2026 (b), measured): on the complete coordinate faces
% (hinf_rd1_hv, stab_rd1_hv) cuADMM needs 35-49x fewer iterations, but the
% reduced programs still have no interior (margin t* ~ 0).  The dual stage did
% give one (hinf_rd1, t* = 1.8e-5), yet its certified solves failed, partly in
% the harness.  Stage 2 needs a dual y, i.e. a full solve, and a dense V fills
% At; no scalable face was found in 2-D.  Closed as a solver-side step; kept
% for reference.
%
% ROWS.  X_k = V_k Y_k V_k', row i -> sum_k <V_k' A_ik V_k, Y_k> + free part.  Rows
% that vanish (b_i = 0 and reduced norm <= drop x original norm) are dropped as
% implied by the face.  A wrong face cannot give a wrong certificate: the
% mapped-back point X = V Y V' (PSD by construction when Y is) is re-checked on
% EVERY original row, dropped ones included.
%
% COST.  kron(V_k,V_k)'*A_k per block: sparse when V_k is a pure selection; dense
% N^2 x N'^2 otherwise, so stage 2 is limited to blocks of up to opts.nmax (100).
%
% OUTPUT Rd: At, b (reduced rows only), K, keeprows, V, Ks, Kr, nzero (stage 1),
% nexp (stage 2), xfull = @(xred_fullvec) original full-vec point.

if nargin < 3, opts = struct(); end
if ~isfield(opts,'sig'),  opts.sig = 0.05; end
if ~isfield(opts,'drop'), opts.drop = 1e-9; end
if ~isfield(opts,'nmax'), opts.nmax = 100; end
if ~isfield(opts,'taux'), opts.taux = 1e-6; end
K = D.K;  if ~isfield(K,'f') || isempty(K.f), K.f = 0; end
Ks = double(K.s(:)');  nb = numel(Ks);  b = full(D.b(:));
F = bl_face(D.At,b,K);
z = [];
if isfield(ref,'y') && ~isempty(ref.y)
    y = ref.y(:);  z1 = -(D.At*y);
    if blockmin(z1,K) < blockmin(-z1,K), z1 = -z1; end        % the sign that is ~PSD
    z = z1;
end
zglob = 0;  off = K.f;                                         % largest |eig(Z)| over all blocks
if ~isempty(z)
    for k = 1:nb
        N = Ks(k);  Zk = reshape(z(off+(1:N^2)),N,N);  off = off + N^2;
        zglob = max(zglob,max(abs(eig((Zk+Zk')/2))));
    end
end
hasx = isfield(ref,'x') && ~isempty(ref.x);  if hasx, xr0 = ref.x(:); end
V = cell(1,nb);  Kr = zeros(1,nb);  nexp = zeros(1,nb);
parts = {D.At(1:K.f,:)};  off = K.f;
for k = 1:nb
    N = Ks(k);  rows = off + (1:N^2);
    kp = find(~F.zero{k});  P = sparse(kp,1:numel(kp),1,N,numel(kp));
    Vk = P;  dense = false;
    if ~isempty(z) && N <= opts.nmax && ~isempty(kp)
        Zk = reshape(z(rows),N,N);  Zk = (Zk+Zk')/2;
        [E,L] = eig(full(P'*Zk*P));  L = diag(L);
        ex = L > opts.sig*zglob;
        if hasx && any(ex)                                     % complementarity: X small along it
            Xk = reshape(xr0(rows),N,N);  Xk = full(P'*((Xk+Xk')/2)*P);
            xl = sum(E.*(Xk*E),1)';  ex = ex & (xl <= opts.taux*max(eig(Xk)));
        end
        nexp(k) = nnz(ex);
        if nexp(k) > 0, Vk = P*E(:,~ex);  dense = true; end
    end
    off = off + N^2;
    if ~dense && numel(kp) == N                                % nothing removed
        V{k} = [];  Kr(k) = N;  parts{end+1} = D.At(rows,:);   %#ok<AGROW>
        continue
    end
    V{k} = Vk;  Kr(k) = size(Vk,2);
    parts{end+1} = sparse(kron(Vk,Vk)'*D.At(rows,:));          %#ok<AGROW>
end
keepblk = Kr > 0;
Ar = vertcat(parts{[true keepblk]});
nr = full(sqrt(sum(Ar.^2,1)))';  n0 = full(sqrt(sum(D.At.^2,1)))';
keeprows = ~(b == 0 & nr <= opts.drop*n0);
Kred = struct('f',K.f,'l',0,'q',[],'s',Kr(keepblk));
Rd = struct('At',Ar(:,keeprows),'b',b(keeprows),'K',Kred,'keeprows',keeprows,'V',{V}, ...
            'Ks',Ks,'Kr',Kr,'nzero',F.nzero,'nexp',nexp,'ndrop',nnz(~keeprows));
Rd.xfull = @(xr) expand(xr,Rd,K);
if isfield(opts,'outmat') && ~isempty(opts.outmat)
    S = struct('At',Rd.At,'b',Rd.b,'c',sparse(size(Rd.At,1),1),'K',Kred,'Ns',Kred.s,'Kf',K.f);
    save(opts.outmat,'-struct','S','-v7.3');
    RR = speye(size(Rd.At,1));  bscl = 1;  Sshape = struct();      %#ok<NASGU>
    save(strrep(opts.outmat,'.mat','_meta.mat'),'RR','bscl','Sshape');
end
end


function m = blockmin(z,K)
m = inf;  off = K.f;
for k = 1:numel(K.s)
    N = K.s(k);  M = reshape(z(off+(1:N^2)),N,N);  off = off + N^2;
    e = eig((M+M')/2);  m = min(m,min(e)/max(max(abs(e)),realmin));
end
end


function x = expand(xr,Rd,K)
% full-vec original point from a full-vec reduced one: X_k = V_k Y_k V_k'
xr = xr(:);  x = zeros(K.f + sum(Rd.Ks.^2),1);  x(1:K.f) = xr(1:K.f);
o = K.f;  orr = K.f;
for k = 1:numel(Rd.Ks)
    N = Rd.Ks(k);  n = Rd.Kr(k);
    if n > 0
        Y = reshape(xr(orr+(1:n^2)),n,n);  orr = orr + n^2;
        if isempty(Rd.V{k}), X = Y; else, X = full(Rd.V{k}*Y*Rd.V{k}'); end
    else
        X = zeros(N);
    end
    x(o+(1:N^2)) = X(:);  o = o + N^2;
end
end
