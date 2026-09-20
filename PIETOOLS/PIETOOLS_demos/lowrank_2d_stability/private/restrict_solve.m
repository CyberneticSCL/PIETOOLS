function [R,q,info] = restrict_solve(prog,H,P,Atf,bf,V)
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
% WHY NOT LEAST SQUARES.  ||b|| ~ 5e-6 on these programs, so X = 0 nearly
% satisfies the equalities; MEASURED on this 2D family (v2d_sweep.m) an lsqr
% point reached rel_eq = 1.7e-11 AND op rel = 1.7e-11 yet had mineig = -4.4e-07,
% i.e. it satisfied the identity while being INDEFINITE -- no certificate.
% Positivity is the whole content, so the restricted problem is solved as an SDP
% and the gate is opcheck_2d's R.ok = (rel < 1e-6) AND psd.
%
% MEMORY.  The 1D version of this code did full(Atf(rows_i,:)); at N = 424 that
% block is 179776 x 3456 = 4.9 GB.  Here V_i' is applied to the SPARSE block on
% both sides instead, so the largest intermediate is r*N*m doubles.
%
% INPUT   prog,H  from build_stab_2d_st2;  P  bm_setup package
%         Atf,bf  ORIGINAL (un-normalised) SeDuMi data from raw_data
%         V       1 x B cell, V{i} is Ns(i) x r_i with orthonormal columns
% OUTPUT  R    opcheck_2d result for the LIFTED solution (R.ok is the verdict)
%         q    lifted Gram vector in NORMALISED-b units (opcheck convention)
%         info .rr .np .rM .t .rel_eq_full .mineig_lift .sedumi_iter

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
for i = 1:B
    ri = rr(i);
    Ab = full(Atr(roff2+(1:ri^2),:));                 % ri^2 x m, ri is tiny
    for a = 1:ri
        for bq = a:ri
            pcol = pcol + 1;
            if a==bq
                M(:,pcol) = Ab((a-1)*ri+a,:)';
            else
                M(:,pcol) = (Ab((bq-1)*ri+a,:) + Ab((a-1)*ri+bq,:))';
            end
        end
    end
    roff2 = roff2 + ri^2;
end

% ---- reduce to independent rows, then the tiny SDP ------------------------
% SeDuMi needs more variables than constraints and here np << m, but
% rank(M) <= np so all but a handful of rows are redundant.  Pivoted QR on M'
% picks an independent set; the verdict is still taken against ALL the original
% rows (and then the operators), so a wrong reduction cannot pass.
bfull = full(bf(:));
[~,Rq,eq_] = qr(full(M'),0);
dq = abs(diag(Rq));
rM = max(nnz(dq > max(dq)*1e-12),1);
keep = sort(eq_(1:rM));
npad = max(0, rM - (Kf + sum(rr.^2)) + 1);            % free padding columns:
Atr_s = [sparse(npad,m); Atr];                        % all-zero, so they cannot
Kr = struct('f',Kf+npad,'l',0,'s',rr);                % change the feasible set
pars.fid = 0;
[xr,~,infr] = sedumi(Atr_s(:,keep),bfull(keep),sparse(size(Atr_s,1),1),Kr,pars);
xr = full(xr);   xr = xr(npad+1:end);
info.t = toc(t0);
info.np = np;  info.rM = rM;  info.sedumi_iter = getf(infr,'iter');

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
