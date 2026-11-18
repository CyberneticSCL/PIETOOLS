% PF_SOL_CHECK  Login-node pre-flight of the Sol package: the path setup,
% the 3-D PIE, and a 1-D box-row program at dp = [0 1] solved with MOSEK
% (seconds). SCRATCH, MMP, 11/15/2026.
S = pf_setup();
pie = heatNd_pie(3,0,{'DD','DN','DN'});   fprintf('3-D PIE ok, lambda1 = %.6f\n',pie.exact.lambda1);
pie1 = heatNd_pie(1,0,{'DD'});   lam1 = pie1.exact.lambda1;
[rows,offs,pd] = pf_terms('box',1,1);
opts = struct('preset','bench','dp',[0 1],'psatz','rows','psatz_rows',rows,'psatz_rows_offset',offs,'psatz_rows_perdir',pd);
[prog,meta] = heatNd_lpi(pie1,0,0.9*lam1,0.1,opts);
v = heatNd_solve(prog,struct('tight',true,'echo',0,'param',struct('MSK_IPAR_NUM_THREADS',2)));
fprintf('@@ CHECK 1-D box o1 dp=[0 1] kappa=0.9: st=%+d relb=%.1e m=%d nx=%d it=%d thr=%d %s\n',v.st,v.rel_b,v.m,v.nx,v.iter,v.mosek_threads,v.why);
assert(v.st==1,'pre-flight: the 1-D solve did not certify');
disp('## CHECK DONE');
