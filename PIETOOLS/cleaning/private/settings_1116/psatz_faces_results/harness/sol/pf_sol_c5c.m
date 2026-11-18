% PF_SOL_C5C  Sol job: the 3-D box set at the next cap (C5c of the
% psatz_faces run, 11/15/2026). Heat DD x DN x DN on [0,1]^3 at kappa =
% 0.985 kappa* (kappa* = 1.5 pi^2), R at the bench caps and Q's plain term
% at the heavy caps (RQ 'Tspan', dp = [0 1]), the box rows {g_1, g_2, g_3}
% at per-direction offset 1 on both R and Q. One certified solve
% (heatNd_solve, tight, MOSEK log on); the independent-row fallback only if
% the verdict is uncertain and under 6 h have elapsed. On the workstation
% the program (m 77,157, 29.6M PSD variables) was built in 295 s and not
% solved (gate 4.5M); the faces at offset 0 certify this kappa at the bench
% caps (2.71M variables) and every box configuration there is Farkas-
% infeasible. Prints '@@ C5C' lines; saves results/c5c_sol.mat.
% SCRATCH, MMP, 11/15/2026.
tstart = tic;
S = pf_setup();
nt = maxNumCompThreads;
hwm = @() str2double(regexp(fileread('/proc/self/status'),'VmHWM:\s*(\d+)','tokens','once'))/2^20;   % peak RSS, GiB
fprintf('@@ C5C start %s node %s threads %d\n',datestr(now,31),getenv('SLURMD_NODENAME'),nt);
pie = heatNd_pie(3,0,{'DD','DN','DN'});    lam1 = pie.exact.lambda1;
k = 0.985*lam1;
[rows,offs,pd] = pf_terms('box',3,1);
opts = struct('preset','bench','dp',[0 1],'psatz','rows','psatz_rows',rows,'psatz_rows_offset',offs,'psatz_rows_perdir',pd);
t0 = tic;   [prog,meta] = heatNd_lpi(pie,0,k,0.1,opts);    tb = toc(t0);
sdp = meta.sdp;
fprintf('@@ C5C box o1 dp=[0 1] build %.0f s: m=%d nx=%d Ks=%s ndec=%d hwm=%.1f GiB\n',tb,sdp.m,sdp.nx,mat2str(sdp.Ks),meta.ndec,hwm());
assert(sdp.m==77157 && sdp.nx==29554587,'C5C: the program differs from the workstation build (m %d, nx %d)',sdp.m,sdp.nx);
out = [S '/results/c5c_sol.mat'];
R = struct('cfg','box','o',1,'pd',pd,'dp',[0 1],'frac',0.985,'kappa',k,'lam1',lam1,'m',sdp.m,'nx',sdp.nx,'Ks',sdp.Ks,'ndec',meta.ndec, ...
           'tb',tb,'tc',NaN,'nnzAt',NaN,'ts',NaN,'v',[],'v2',[],'keep',[],'node',getenv('SLURMD_NODENAME'),'threads',nt,'hwm',[]);
t0 = tic;   D = heatNd_sdp(prog);   R.tc = toc(t0);  clear prog
R.nnzAt = nnz(D.At);    R.hwm(end+1) = hwm();
fprintf('@@ C5C sdp %.0f s: nnzAt=%d hwm=%.1f GiB\n',R.tc,R.nnzAt,R.hwm(end));
save(out,'R');
sopts = struct('tight',true,'echo',3,'param',struct('MSK_IPAR_NUM_THREADS',nt,'MSK_IPAR_LOG_INTPNT',1,'MSK_DPAR_OPTIMIZER_MAX_TIME',30600));
t0 = tic;   v = heatNd_solve(D,sopts);    R.ts = toc(t0);   R.v = rmfield(v,'x');  R.hwm(end+1) = hwm();  save(out,'R');
fprintf('@@ C5C box o1 dp=[0 1] kappa=%.4f (0.985) st=%+d relb=%.1e psd=%.1e viol=%.1e m=%d nx=%d it=%d tmosek=%.0f ts=%.0f thr=%d hwm=%.1f GiB %s | %s / %s\n', ...
        k,v.st,v.rel_b,v.psd_relmin,v.cert_viol,v.m,v.nx,v.iter,v.t_mosek,R.ts,v.mosek_threads,R.hwm(end),v.why,v.prosta,v.solsta);
if v.st==0 && toc(tstart) < 6*3600
    t0 = tic;   [keep,li] = heatNd_lindep(D.At);  tl = toc(t0);
    fprintf('@@ C5C lindep %.0f s: rank %d of %d hwm=%.1f GiB\n',tl,li.rank,li.m,hwm());
    sopts.rows = 'keep';    sopts.keep = keep;
    sopts.param.MSK_DPAR_OPTIMIZER_MAX_TIME = max(600, 9*3600 - toc(tstart));
    t0 = tic;   v2 = heatNd_solve(D,sopts);  ts2 = toc(t0);  R.v2 = rmfield(v2,'x');   R.keep = keep;   R.hwm(end+1) = hwm();  save(out,'R');
    fprintf('@@ C5C box o1 dp=[0 1] keep kappa=%.4f (0.985) st=%+d relb=%.1e psd=%.1e viol=%.1e nrows=%d it=%d ts=%.0f hwm=%.1f GiB %s | %s / %s\n', ...
            k,v2.st,v2.rel_b,v2.psd_relmin,v2.cert_viol,v2.nrows,v2.iter,ts2,R.hwm(end),v2.why,v2.prosta,v2.solsta);
end
fprintf('## C5C DONE %.0f s\n',toc(tstart));
