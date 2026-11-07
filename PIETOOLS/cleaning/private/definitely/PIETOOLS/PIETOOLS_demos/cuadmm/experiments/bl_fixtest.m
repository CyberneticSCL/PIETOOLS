% bl_fixtest.m -- does pinning gamma at the SDP level actually pose the
% feasibility question, and what separates feasible from infeasible?
%
% hinf_rd1's objective form returned gam* = 0.182626 under Mosek. If the
% transformation is right, the pinned program must be solvable for gam >= gam*
% and not for gam < gam*, and the transition must be sharp enough to bisect on.
% Swept well to either side so the cliff is located rather than assumed, and
% every quantity that a decision rule might use is printed side by side --
% rel_b, Mosek's status words, the PSD margin and ||x|| -- because in this
% campaign numerr has already passed a program with rel_b = 1.000 and feasratio
% has read -0.06 at a comfortably feasible point.
cuadmm_path;
HERE = fileparts(mfilename('fullpath'));
DMP  = fullfile(cuadmm_outdir(),'baseline','dumps');
ID   = 'hinf_rd1';
GSTAR = 0.182626;
FR = [0.50 0.80 0.95 0.99 0.999 1.0 1.001 1.01 1.05 1.20 2.00];

fprintf('FX frac|gam|m|rel_b|psd_min|normx|prosta|solsta|t_s\n');
for f = FR
    g = f*GSTAR;
    D = bl_fix(fullfile(DMP,[ID '.mat']),g);
    prob = Sedumi2Mosek(D.At',full(D.b),D.c,D.K);
    t0 = tic; [~,res] = mosekopt('minimize info echo(0)',prob); ts = toc(t0);
    x  = MosekSol2SedumiSol(D.K,res);  x = x(:);
    rb = norm(full(D.At'*x - D.b))/max(norm(full(D.b)),eps);
    off = D.K.f; pmin = inf; pmax = -inf;
    for k = 1:numel(D.Ks)
        N = double(D.Ks(k));
        Xk = reshape(x(off+(1:N^2)),N,N); ev = eig((Xk+Xk')/2);
        pmin = min(pmin,min(ev)); pmax = max(pmax,max(ev)); off = off + N^2;
    end
    fprintf('FX %.3f|%.6f|%d|%.4e|%+.3e|%.4e|%s|%s|%.3f\n', ...
        f,g,D.m,rb,pmin/max(pmax,eps),norm(x), ...
        res.sol.itr.prosta,res.sol.itr.solsta,ts);
end
fprintf('FXDONE\n');
