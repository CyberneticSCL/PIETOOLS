% PF_RUN_C5C  3-D: the box set with Q's plain term one degree up (RQ
% 'Tspan', dp = [0 1]: R at the bench caps, Q at the heavy caps), per-
% direction offsets, at kappa = 0.985 kappa*, where faces0 at dp 0
% certifies (2.71M variables) and box0 at dp 0 is Farkas-infeasible. Size-
% gated: the program is built and its size printed; it is solved only
% below 4.5M PSD variables. SCRATCH, MMP, 11/15/2026.
S = pf_setup();
pie = heatNd_pie(3,0,{'DD','DN','DN'});    lam1 = pie.exact.lambda1;
k = 0.985*lam1;
[rows,offs,pd] = pf_terms('box',3,1);
opts = struct('preset','bench','dp',[0 1],'psatz','rows','psatz_rows',rows,'psatz_rows_offset',offs,'psatz_rows_perdir',pd);
t0 = tic;   [prog,meta] = heatNd_lpi(pie,0,k,0.1,opts);    tb = toc(t0);
sdp = meta.sdp;
fprintf('@@ C5C box o1 dp=[0 1] build %.0f s: m=%d nx=%d Ks=%s ndec=%d\n',tb,sdp.m,sdp.nx,mat2str(sdp.Ks),meta.ndec);
if sdp.nx <= 4.5e6
    t0 = tic;   v = heatNd_solve(prog,struct('tight',true));    ts = toc(t0);
    if v.st==0
        D = heatNd_sdp(prog);   D.keep = heatNd_lindep(D.At);
%       t0 = tic;   v = heatNd_solve(D,struct('tight',true,'rows','keep'));  ts = ts + toc(t0); % MMP, 11/15/2026 (was: no OPTS.keep, would error; D31)
        t0 = tic;   v = heatNd_solve(D,struct('tight',true,'rows','keep','keep',D.keep));  ts = ts + toc(t0); % MMP, 11/15/2026
    end
    mem = NaN;  if isfield(v,'mem_peak'),   mem = v.mem_peak;   end
    fprintf('@@ C5C box o1 dp=[0 1] kappa=%.4f (0.985) st=%+d relb=%.1e psd=%.1e viol=%.1e m=%d nx=%d ts=%.0f memMB=%.0f %s\n', ...
            k,v.st,v.rel_b,v.psd_relmin,v.cert_viol,v.m,v.nx,ts,mem,v.why);
else
    fprintf('@@ C5C box o1 dp=[0 1]: %d PSD variables, above the 4.5M gate; not solved\n',sdp.nx);
end
disp('## C5C DONE');
