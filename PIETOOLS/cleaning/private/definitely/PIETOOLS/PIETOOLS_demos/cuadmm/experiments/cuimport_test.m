% cuimport_test.m -- validate the import path BEFORE any number rests on it.
%
% Three things are established here, in order of what they rule out:
%
% 1. FORWARD check of T on known-good data.  Mosek's cone vector x_m solves the
%    dumped system.  T*x_m must then satisfy cuADMM's system, because T'*T is the
%    identity on the symmetric subspace.  This exercises T in the EXPORT
%    direction on a vector T did not produce, so it is not self-referential.
%
% 2. IMPORT of cuADMM's own X, scored on the original data.
%
% 3. NEGATIVE CONTROL.  The import is re-run with a deliberately wrong
%    off-diagonal factor.  The claim being tested is that the row residual is
%    BLIND to an inverse-pair error while the PSD margin is not.  If the
%    corrupted import also shows a healthy psd_min, then the PSD check has no
%    power and check 2 proves nothing -- exactly the degenerate-test failure that
%    cost a wrong conclusion earlier today.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
D  = fullfile(cuadmm_outdir(),'cu_fish','fish_nobnd');
MF = fullfile(cuadmm_outdir(),'cu_fish','fish_nobnd.mat');
S  = load(MF);

% ---- 1. forward check, needs a Mosek point on the SAME dumped system
t0 = tic; pars.fid = 0;
[xm,ym,im] = sedumi(S.At,full(S.b),S.c,S.K,pars);     % SeDuMi: no Mosek bridge for raw cone data
tm = toc(t0);
fwd_rel = norm(full(S.At'*xm - S.b(:)))/norm(full(S.b(:)));
K = S.K; Ks = double(K.s(:)'); Kf = K.f;
vec_len = Kf + sum(Ks.*(Ks+1)/2);
T = build_T(Kf,Ks,vec_len,size(S.At,1),sqrt(2)/2);
cu_At = T*S.At;  sx = T*xm;
fwd_cu = norm(full(cu_At'*sx - S.b(:)))/norm(full(S.b(:)));
fprintf('CI fwd|sedumi_rel_b=%.4e|cuADMM_system_rel_b=%.4e|t=%.1f|numerr=%d\n', ...
        fwd_rel,fwd_cu,tm,im.numerr);

% ---- 2. import cuADMM's answer
r = cuimport(D,MF);
fprintf('CI imp|rel_b=%.4e|psd_min=%+.4e|psd_relmin=%+.4e|max_asym=%.2e|dual_min=%+.4e|obj=%.6e|normx=%.4e\n', ...
        r.rel_b_norm,r.psd_min,r.psd_relmin,r.max_asym,r.dual_min,r.obj_norm,r.normx);
fprintf('CI blk|%s\n',strjoin(arrayfun(@(v)sprintf('%+.2e',v),r.blk_mineig,'uni',0),' '));

% ---- 3. negative control: wrong off-diagonal factor on the IMPORT side only
xs = fscanf(fopen(fullfile(D,'X_opt.txt'),'r'),'%f');
for fac = [sqrt(2)/2, 0.5, 1.0]
    Tw = build_T(Kf,Ks,vec_len,size(S.At,1),fac);
    xw = Tw'*xs(:);
    rbw = norm(full(S.At'*xw - S.b(:)))/norm(full(S.b(:)));
    off = Kf; pm = inf; px = -inf;
    for k = 1:numel(Ks)
        N = Ks(k); Xk = reshape(xw(off+(1:N^2)),N,N);
        ev = eig((Xk+Xk')/2); pm = min(pm,min(ev)); px = max(px,max(ev));
        off = off + N^2;
    end
    fprintf('CI neg|factor=%.4f|rel_b=%.4e|psd_min=%+.4e|psd_relmin=%+.4e\n', ...
            fac,rbw,pm,pm/max(px,eps));
end
fprintf('CIDONE\n');

function T = build_T(Kf,Ks,vec_len,nvar,r2)
I={};J={};V={};
if Kf>0, p=(1:Kf)'; I{end+1}=p; J{end+1}=p; V{end+1}=ones(Kf,1); end
bs=Kf; bv=Kf;
for k=1:numel(Ks)
    N=Ks(k); idx=find(triu(true(N))); [ii,jj]=ind2sub([N N],idx);
    p=(1:numel(idx))'; dg=(ii==jj); od=~dg;
    I{end+1}=[bs+p; bs+p(od)];                                       %#ok<AGROW>
    J{end+1}=[bv+idx; bv+sub2ind([N N],jj(od),ii(od))];              %#ok<AGROW>
    V{end+1}=[dg + od*r2; repmat(r2,nnz(od),1)];                     %#ok<AGROW>
    bs=bs+N*(N+1)/2; bv=bv+N^2;
end
T=sparse(cat(1,I{:}),cat(1,J{:}),cat(1,V{:}),vec_len,nvar);
end
