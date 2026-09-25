% cu_hinf_gam.m -- convert cuADMM's objective to gamma and compare.
% b was normalised by bscl, so the argmin scaled by 1/bscl and gamma scales
% BACK: gamma = pobj * bscl.  Cross-checked against the surviving X_opt.txt
% imported through the svec inverse and RRx (only the LAST tolerance's
% solution survives per directory -- the stale-solution hazard).
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
ROOT = fullfile(cuadmm_outdir(),'cu_hinf');
% pobj as printed by cuADMM
R = { 'light',1e-4,7.38625754e-02,31859 ,39.54
      'light',1e-6,7.38792316e-02,200000,252.81
      'light',1e-8,7.38792316e-02,200000,450.98
      'heavy',1e-4,7.37608386e-02,3074  ,4.89
      'heavy',1e-6,7.37076882e-02,78105 ,127.29
      'heavy',1e-8,7.37075706e-02,200000,325.46 };
fprintf('HG set|tol|iters|t|gamma_cuadmm|gamma_ref|rel_err\n');
for k = 1:size(R,1)
    M = load(fullfile(ROOT,['hinf_' R{k,1} '_meta.mat']));
    g = R{k,3}*M.bscl;
    fprintf('HG %-5s|%.0e|%d|%.1f|%.9f|%.9f|%.2e\n', ...
        R{k,1},R{k,2},R{k,4},R{k,5},g,M.g_ref,abs(g-M.g_ref)/abs(M.g_ref));
end

fprintf('HG --- cross-check: import the surviving X_opt through svec-inverse + RRx\n');
for pre = {'light','heavy'}
    M = load(fullfile(ROOT,['hinf_' pre{1} '_meta.mat']));
    f = fullfile(ROOT,['hinf_' pre{1}],'X_opt.txt');
    if ~exist(f,'file'), fprintf('HG %s no X_opt\n',pre{1}); continue; end
    v = load(f);  v = v(:);
    x = unsvec(v, M.S.Kf, M.S.Ks);
    RRx = (M.RR*x)*M.bscl;
    fprintf('HG %-5s|import gamma=%.9f | ref=%.9f | rel=%.2e\n', ...
        pre{1}, RRx(M.obj_idx), M.g_ref, abs(RRx(M.obj_idx)-M.g_ref)/abs(M.g_ref));
end
fprintf('HGDONE\n');

function x = unsvec(v,Kf,Ks)
% inverse of the SDPT3 svec used by dump2cuadmm: upper triangle column by
% column, diagonal unscaled, off-diagonal scaled by sqrt(2).
x = v(1:Kf);  p = Kf;
for k = 1:numel(Ks)
    N = Ks(k);  X = zeros(N);
    for j = 1:N
        for i = 1:j
            p = p+1;
            if i==j, X(i,j) = v(p); else, X(i,j) = v(p)/sqrt(2); X(j,i) = X(i,j); end
        end
    end
    x = [x; X(:)];
end
end
