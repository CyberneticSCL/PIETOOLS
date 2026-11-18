% PF_RUN_C5B  Phase 4b again: the 3-D heat benchmark at the bench degrees
% with the weighted terms at offset 0 (the record's configuration: faces
% at the full degree certified [14.60, 14.70) of 14.804): faces0, box0 and
% product0 at kappa = 0.985 kappa*, then box0 at 0.992 kappa*. One
% certified solve each; sizes, times, peak memory. SCRATCH, MMP, 11/15/2026.
S = pf_setup();
pie = heatNd_pie(3,0,{'DD','DN','DN'});    lam1 = pie.exact.lambda1;
out = sprintf('%s/results/c5b_3d.mat',S);
L = {'faces0',0.985; 'box0',0.985; 'product0',0.985; 'box0',0.992};
R = struct('cfg',{},'frac',{},'st',{},'why',{},'rel_b',{},'psd',{},'viol',{},'m',{},'nx',{},'Ks',{},'ndec',{},'tb',{},'ts',{},'mem',{},'rows',{},'err',{});
for i = 1:size(L,1)
    s = L{i,1};     k = L{i,2}*lam1;
    [rows,offs,pd] = pf_terms(s,3,[]);
    r = struct('cfg',s,'frac',L{i,2},'st',NaN,'why','','rel_b',NaN,'psd',NaN,'viol',NaN,'m',NaN,'nx',NaN,'Ks',[],'ndec',NaN,'tb',NaN,'ts',NaN,'mem',NaN,'rows','all','err','');
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
            s,k,L{i,2},r.st,r.rel_b,r.psd,r.viol,r.m,r.nx,mat2str(r.Ks),r.ndec,r.tb,r.ts,r.mem,r.rows,r.why,r.err);
    save(out,'R');
end
disp('## C5B DONE');
