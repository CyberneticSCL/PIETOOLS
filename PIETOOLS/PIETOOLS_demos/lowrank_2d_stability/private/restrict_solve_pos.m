function [R,q,info] = restrict_solve_pos(prog,Ptgt,Peop,P,Atf,bf,V)
% restrict_solve_pos -- restrict_solve.m for BARE POSITIVITY programs: the only
% change is the acceptance gate, opcheck_pos2d (Ptgt vs Peop(q)) instead of
% opcheck_2d's stability identity.  The restriction/lift mechanics, the
% sparse-side application of V', the (:) find-orientation trap at r = 1, the
% pivoted-QR row reduction and the free-padding for SeDuMi's variable count are
% copied verbatim from the validated restrict_solve.m -- see its header for
% the soundness (one-way) and memory arguments.
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
    [ii,jj,vv] = find(sparse(Atil));
    Ri{end+1}=ii(:)+roff; Rj{end+1}=jj(:); Rv{end+1}=vv(:); %#ok<AGROW>
    off = off + N^2;   roff = roff + ri^2;
end
Atr = sparse(cat(1,Ri{:}),cat(1,Rj{:}),cat(1,Rv{:}),Kf+sum(rr.^2),m);

% ---- M: the restricted map in SYMMETRIC parameters ------------------------
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
bfull = full(bf(:));
[~,Rq,eq_] = qr(full(M'),0);
dq = abs(diag(Rq));
rM = max(nnz(dq > max(dq)*1e-12),1);
keep = sort(eq_(1:rM));
npad = max(0, rM - (Kf + sum(rr.^2)) + 1);
Atr_s = [sparse(npad,m); Atr];
Kr = struct('f',Kf+npad,'l',0,'s',rr);
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

% opcheck wants NORMALISED-b units (trap 2: omitting nb0 reports rel = 1)
q = zeros(P.Ntot,1);   off = 0;
for i = 1:B
    N = Ns(i);
    q(P.rows{i}) = xfull(Kf+off+(1:N^2))/P.nb0;   off = off+N^2;
end
if Kf>0, q(1:Kf) = xfull(1:Kf)/P.nb0; end
R = opcheck_pos2d(prog,Ptgt,Peop,P,q);
end

function v = getf(s,f)
if isstruct(s) && isfield(s,f), v = s.(f); else, v = -1; end
end
