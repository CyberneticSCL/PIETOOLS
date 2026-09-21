function [R,q,info] = restrict_solve(prog,H,P,Atf,bf,V,use_sdp)             % CC, 09/20/2026
% restrict_solve -- restrict the SDP to the face X_i = V_i S_i V_i', solve the
% tiny problem, lift, and VERIFY AT THE OPERATOR LEVEL.  This is the acceptance
% test for a claimed rank; nothing else in this study may accept one.
%
% SOUNDNESS, one way.  Any S_i >= 0 gives X_i = V_i S_i V_i' >= 0, so the
% restriction can only LOSE feasibility, never manufacture a false certificate.
% A rank accepted here is therefore a genuine upper bound on the minimum rank,
% and a rank rejected here is only "not reachable on THIS face" -- stated as
% such, never as a proven lower bound.
%
% WHY NOT LEAST SQUARES ON THE FULL PROGRAM.  ||b|| ~ 5e-6 on these programs,
% so X = 0 nearly satisfies the equalities; MEASURED on this 2D family
% (v2d_sweep.m) an lsqr point reached rel_eq = 1.7e-11 AND op rel = 1.7e-11 yet
% had mineig = -4.4e-07, i.e. it satisfied the identity while being INDEFINITE
% -- no certificate.  Positivity is the whole content, so the gate is
% opcheck_2d's R.ok = (rel < 1e-6) AND psd, never an equality residual.
% That trap is a property of the FULL program, which is massively
% underdetermined: min-norm is free to pick an indefinite point.  It does not
% carry over to the RESTRICTED system -- see the next paragraph.
%
% HOW THE RESTRICTED SYSTEM IS SOLVED (CC, 09/20/2026).  On the face, the kept
% rows M(keep,:)*s = b(keep) have FULL COLUMN RANK at every working rank
% measured here (np - rM = 0 for r = 2..30, sigma_min 0.11-0.99 at r = 1..15),
% so the restricted problem is a DETERMINED linear solve, not an optimisation:
% exactly one s meets those rows and any solver must return it.  Handing it to
% an interior-point solver therefore bought nothing and could fail -- measured,
% SeDuMi reported pinf and returned x = 0 exactly on faces that CONTAIN a
% certificate (a rank-5 face padded to r = 6/7/10/15/20/30, and a lifted face
% at cond(M) = 3.5e16), i.e. FALSE NEGATIVES on the one verdict this file
% exists to produce.  s is therefore taken from the backslash solve, and S is
% then PROJECTED ONTO THE PSD CONE: the determined solve is sign-free, so its
% unique point carries small negative eigenvalues (measured -2.6e-08 to
% -5.2e-06 relative) in directions the face barely uses, and without the
% projection the operator residual is repaired (1 -> 1e-07) while the PSD half
% of the gate still rejects.  SeDuMi is not "wrong" on those faces: no PSD
% point satisfies the kept rows EXACTLY.  Acceptance never asked for that --
% it asks for a PSD Gram whose operator residual is below the gate.
% None of this weakens acceptance.  The solve asserts nothing; the gate below
% -- opcheck's 36-cell operator residual AND per-block PSD, unchanged -- was
% always the only thing allowed to accept a rank, and it re-measures the
% operator residual AT the projected point, so a projection large enough to
% matter shows up as a failed gate (measured: a random face whose solution
% clips from mineig/maxeig = -0.94 still reports op rel 0.995, rejected).
% Where rM < np the backslash returns a basic (not min-norm) solution and the
% gate is likewise the only verdict.
% use_sdp = true keeps the SeDuMi path callable as the oracle.
%
% MEMORY.  The 1D version of this code did full(Atf(rows_i,:)); at N = 424 that
% block is 179776 x 3456 = 4.9 GB.  Here V_i' is applied to the SPARSE block on
% both sides instead, so the largest intermediate is r*N*m doubles.
%
% INPUT   prog,H  from build_stab_2d_st2;  P  bm_setup package
%         Atf,bf  ORIGINAL (un-normalised) SeDuMi data from raw_data
%         V       1 x B cell, V{i} is Ns(i) x r_i with orthonormal columns
%         use_sdp optional, default false: true solves the restricted system
%                 with SeDuMi as before -- the oracle path, and the only
%                 reason to set it
% OUTPUT  R    opcheck_2d result for the LIFTED solution (R.ok is the verdict)
%         q    lifted Gram vector in NORMALISED-b units (opcheck convention)
%         info .rr .np .rM .t .rel_eq_full .mineig_lift .sedumi_iter .method
%              .clip  largest relative negative eigenvalue projected away
%                     (NaN on the SeDuMi path).  A large .clip on an ACCEPTED
%                     point is worth reading: the gate passed it, but the face
%                     was carrying that much indefiniteness before projection.
%
% CC, 09/20/2026: solve the restricted system by backslash on the kept rows
%          and project S onto the PSD cone, instead of calling SeDuMi: the
%          system is determined (np == rM measured at every working rank) and
%          the SDP call returned x = 0 on faces that contain a certificate.
%          MEASURED, both arms in one session on identical faces: no change on
%          8 banked certificates (op rel agrees to every printed digit), and
%          7 faces that SeDuMi rejected while containing a certificate are now
%          accepted; 9 negative controls stay rejected.  The gate is untouched.
%          This supersedes the header's "WHY NOT LEAST SQUARES" paragraph AS
%          IT APPLIES TO THE RESTRICTED SOLVE ONLY; that paragraph's
%          measurement (an lsqr point on the FULL program, indefinite at op
%          rel 1.7e-11) still stands and is still the reason the gate, not the
%          equality residual, decides.

if nargin<7 || isempty(use_sdp), use_sdp = false; end                       % CC, 09/20/2026
Ns = P.Ns;   Kf = P.Kf;   m = P.m;   B = numel(Ns);
rr = zeros(1,B);
for i = 1:B, rr(i) = size(V{i},2); end
info.rr = rr;
t0 = tic;

% ---- restrict: Atil(:,k) = vec( V_i' A_k^(i) V_i ) ------------------------
Ri = {};  Rj = {};  Rv = {};  off = 0;  roff = Kf;
if Kf > 0
    [ii,jj,vv] = find(Atf(1:Kf,:));
    Ri{end+1}=ii(:); Rj{end+1}=jj(:); Rv{end+1}=vv(:);
end
for i = 1:B
    N = Ns(i);  ri = rr(i);  Vi = V{i};
    Ab = Atf(Kf+off+(1:N^2),:);                       % N^2 x m, SPARSE, kept so
    M1 = Vi' * reshape(Ab,N,N*m);                     % r x (N*m)
    T  = reshape(M1,ri,N,m);
    M2 = Vi' * reshape(permute(T,[2 1 3]),N,ri*m);    % r x (r*m)
    Atil = reshape(permute(reshape(M2,ri,ri,m),[2 1 3]),ri^2,m);
    % (:) is REQUIRED: at r = 1 Atil is a ROW vector and MATLAB's find returns
    % row-oriented output for row input, which breaks the cat(1,...) below.
    % r = 1 is exactly the interesting case here (trap 3).
    [ii,jj,vv] = find(sparse(Atil));
    Ri{end+1}=ii(:)+roff; Rj{end+1}=jj(:); Rv{end+1}=vv(:); %#ok<AGROW>
    off = off + N^2;   roff = roff + ri^2;
end
Atr = sparse(cat(1,Ri{:}),cat(1,Rj{:}),cat(1,Rv{:}),Kf+sum(rr.^2),m);

% ---- M: the restricted map in SYMMETRIC parameters ------------------------
% <T,S> = sum_i T_ii S_ii + sum_{i<j} (T_ij + T_ji) S_ij
np = sum(rr.*(rr+1)/2);
M = zeros(m,np);  pcol = 0;  roff2 = Kf;
% pidx{p}: the entries of the xr layout that carry s(p), so the symmetric
% solution can be written back without a second nested loop.  Atr's rows ARE
% the xr layout ([Kf free; vec(S_1); ...]), and roff2 is that offset already.
pidx = cell(1,np);   % tiny: np = sum r(r+1)/2, never a q-scaled dimension  % CC, 09/20/2026
for i = 1:B
    ri = rr(i);
    Ab = full(Atr(roff2+(1:ri^2),:));                 % ri^2 x m, ri is tiny
    for a = 1:ri
        for bq = a:ri
            pcol = pcol + 1;
            if a==bq
                M(:,pcol) = Ab((a-1)*ri+a,:)';
                pidx{pcol} = roff2+(a-1)*ri+a;                              % CC, 09/20/2026
            else
                M(:,pcol) = (Ab((bq-1)*ri+a,:) + Ab((a-1)*ri+bq,:))';
                pidx{pcol} = roff2+[(bq-1)*ri+a, (a-1)*ri+bq];              % CC, 09/20/2026
            end
        end
    end
    roff2 = roff2 + ri^2;
end

% ---- CC, 09/20/2026, START: backslash + PSD projection instead of SeDuMi --
% ---- reduce to independent rows, then the tiny DETERMINED solve -----------
% rank(M) <= np so all but a handful of rows are redundant.  Pivoted QR on M'
% picks an independent set; the verdict is still taken against ALL the original
% rows (and then the operators), so a wrong reduction cannot pass.
% (was: "...then the tiny SDP", opening "SeDuMi needs more variables than
%  constraints and here np << m, but rank(M) <= np ..." -- that requirement
%  left with the SeDuMi call.  The selection itself is unchanged.)
bfull = full(bf(:));
[~,Rq,eq_] = qr(full(M'),0);
dq = abs(diag(Rq));
rM = max(nnz(dq > max(dq)*1e-12),1);
keep = sort(eq_(1:rM));
if use_sdp                             % oracle: the pre-09/20 SeDuMi path  % CC, 09/20/2026
    npad = max(0, rM - (Kf + sum(rr.^2)) + 1);       % free padding columns:
    Atr_s = [sparse(npad,m); Atr];                   % all-zero, so they cannot
    Kr = struct('f',Kf+npad,'l',0,'s',rr);           % change the feasible set
    pars.fid = 0;
    [xr,~,infr] = sedumi(Atr_s(:,keep),bfull(keep),sparse(size(Atr_s,1),1),Kr,pars);
    xr = full(xr);   xr = xr(npad+1:end);
    info.sedumi_iter = getf(infr,'iter');
    info.clip = NaN;                     % the cone was the solver's job    % CC, 09/20/2026
else                                                                        % CC, 09/20/2026
    % free variables ride along in the same solve; Kf = 0 in every program  % CC, 09/20/2026
    % this package builds, and M covers the PSD blocks only                 % CC, 09/20/2026
    s = [full(Atr(1:Kf,keep))', M(keep,:)] \ bfull(keep);                   % CC, 09/20/2026
    xr = zeros(Kf+sum(rr.^2),1);   xr(1:Kf) = s(1:Kf);                      % CC, 09/20/2026
    for p = 1:np, xr(pidx{p}) = s(Kf+p); end         % S_ab = S_ba = s(p)   % CC, 09/20/2026
    % project S onto the PSD cone (header): the sign-free solve leaves      % CC, 09/20/2026
    % small negative eigenvalues in directions the face barely uses.  ri is % CC, 09/20/2026
    % tiny, so B eigendecompositions cost nothing at any decision count.    % CC, 09/20/2026
    roff3 = Kf;   info.clip = 0;                                            % CC, 09/20/2026
    for i = 1:B                                                             % CC, 09/20/2026
        ri = rr(i);                                                         % CC, 09/20/2026
        Si = reshape(xr(roff3+(1:ri^2)),ri,ri);   Si = (Si+Si')/2;          % CC, 09/20/2026
        [Wc,Dc] = eig(Si);   dc = diag(Dc);                                 % CC, 09/20/2026
        if min(dc) < 0                                                      % CC, 09/20/2026
            info.clip = max(info.clip,-min(dc)/max(max(dc),realmin));       % CC, 09/20/2026
            Si = Wc*diag(max(dc,0))*Wc';   xr(roff3+(1:ri^2)) = Si(:);      % CC, 09/20/2026
        end                                                                 % CC, 09/20/2026
        roff3 = roff3 + ri^2;                                               % CC, 09/20/2026
    end                                                                     % CC, 09/20/2026
    info.sedumi_iter = -1;                           % no SeDuMi call       % CC, 09/20/2026
end                                                                         % CC, 09/20/2026
if use_sdp, info.method = 'sdp'; else, info.method = 'ls'; end              % CC, 09/20/2026
% ---- CC, 09/20/2026, END --------------------------------------------------
info.t = toc(t0);
info.np = np;  info.rM = rM;
%info.sedumi_iter = getf(infr,'iter');   % now set in the branch above      % CC, 09/20/2026 (was)

% ---- lift and verify against the FULL ORIGINAL program --------------------
xfull = zeros(Kf+sum(Ns.^2),1);
if Kf>0, xfull(1:Kf) = xr(1:Kf); end
off = 0;  roff = Kf;  mn = zeros(1,B);  mx = zeros(1,B);
for i = 1:B
    N = Ns(i);  ri = rr(i);
    Si = reshape(xr(roff+(1:ri^2)),ri,ri);  Si = (Si+Si')/2;
    Xi = V{i}*Si*V{i}';   Xi = (Xi+Xi')/2;
    xfull(Kf+off+(1:N^2)) = Xi(:);
    ev = eig(Si);  mn(i) = min(ev);  mx(i) = max(ev);
    off = off + N^2;   roff = roff + ri^2;
end
info.rel_eq_full = norm(Atf'*xfull - bfull)/max(norm(bfull),realmin);
info.mineig_lift = min(mn./max(mx,realmin));

% opcheck_2d wants NORMALISED-b units -- it multiplies by P.nb0 internally
% (trap 2).  The restricted solve used the ORIGINAL bf, so divide.  Omitting
% this scales the Gram by 5e-6, leaving Pop ~ eppos*I and reporting rel = 1 for
% a perfectly good solution.
q = zeros(P.Ntot,1);   off = 0;
for i = 1:B
    N = Ns(i);
    q(P.rows{i}) = xfull(Kf+off+(1:N^2))/P.nb0;   off = off+N^2;
end
if Kf>0, q(1:Kf) = xfull(1:Kf)/P.nb0; end
R = opcheck_2d(prog,H,P,q);
end

function v = getf(s,f)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = -1; end
end
