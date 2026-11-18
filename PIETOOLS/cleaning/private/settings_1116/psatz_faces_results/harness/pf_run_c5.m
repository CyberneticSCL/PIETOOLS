% PF_RUN_C5  Phase 4b: the 3-D heat benchmark (DD x DN x DN, kappa* =
% 1.5 pi^2) at the bench degrees (RQ 'Tspan', dp 0, d = 0, ep = 0.1), one
% certified solve per term set at kappa = 0.985 kappa* (inside the bench
% bracket [14.60, 14.70) of the record): box (4 terms), product (8), faces
% (7), each with per-direction offsets; sizes, times and peak memory. An
% uncertain verdict is re-solved on the independent row set. One process
% at a time (about 20 GB). SCRATCH, MMP, 11/15/2026.
S = pf_setup();
pie = heatNd_pie(3,0,{'DD','DN','DN'});    lam1 = pie.exact.lambda1;
k = 0.985*lam1;
out = sprintf('%s/results/c5_3d.mat',S);
R = struct('cfg',{},'st',{},'why',{},'rel_b',{},'psd',{},'viol',{},'m',{},'nx',{},'Ks',{},'ndec',{},'tb',{},'ts',{},'mem',{},'rows',{},'err',{});
for s = {'box','product','faces'}
    [rows,offs,pd] = pf_terms(s{1},3,1);
    r = struct('cfg',s{1},'st',NaN,'why','','rel_b',NaN,'psd',NaN,'viol',NaN,'m',NaN,'nx',NaN,'Ks',[],'ndec',NaN,'tb',NaN,'ts',NaN,'mem',NaN,'rows','all','err','');
    try
        opts = struct('preset','bench','psatz','rows','psatz_rows',rows,'psatz_rows_offset',offs,'psatz_rows_perdir',pd);
        t0 = tic;   [prog,meta] = heatNd_lpi(pie,0,k,0.1,opts);    r.tb = toc(t0);
        r.ndec = meta.ndec;
        t0 = tic;   v = heatNd_solve(prog,struct('tight',true));    r.ts = toc(t0);
        if v.st==0
            D = heatNd_sdp(prog);   D.keep = heatNd_lindep(D.At);
%           t0 = tic;   v = heatNd_solve(D,struct('tight',true,'rows','keep'));  r.ts = r.ts + toc(t0);   r.rows = 'keep'; % MMP, 11/15/2026 (was: no OPTS.keep, would error; D31)
            t0 = tic;   v = heatNd_solve(D,struct('tight',true,'rows','keep','keep',D.keep));  r.ts = r.ts + toc(t0);   r.rows = 'keep'; % MMP, 11/15/2026
        end
        r.st = v.st;    r.why = v.why;  r.rel_b = v.rel_b;  r.psd = v.psd_relmin;  r.viol = v.cert_viol;
        r.m = v.m;      r.nx = v.nx;    r.Ks = v.Ks;
        if isfield(v,'mem_peak'),   r.mem = v.mem_peak;     end
    catch ME
        r.err = strrep(ME.message,newline,' ');
    end
    R(end+1) = r;                                                           %#ok<AGROW>
    fprintf('@@ C5 cfg=%s kappa=%.4f (%.3f) st=%+d relb=%.1e psd=%.1e viol=%.1e m=%d nx=%d Ks=%s ndec=%d tb=%.0f ts=%.0f memMB=%.0f rows=%s %s %s\n', ...
            s{1},k,k/lam1,r.st,r.rel_b,r.psd,r.viol,r.m,r.nx,mat2str(r.Ks),r.ndec,r.tb,r.ts,r.mem,r.rows,r.why,r.err);
    save(out,'R');
end
disp('## C5 DONE');
