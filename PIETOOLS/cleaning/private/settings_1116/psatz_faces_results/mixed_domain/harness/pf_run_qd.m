% PF_RUN_QD  The direct PIE2PDE form (PF_QDIRECT) on S3 (0.9 r*, control),
% S17 (0.5 r2*) and S18 (0.5 tau*): Q at degree 1 (every cap 1: mult, int,
% out, in) and at degree 2; dw 0 and 1; eppos 1e-4, epneg 0. SCRATCH, MMP,
% 11/16/2026.
tau43 = 1.034726;
qdegs = {struct('mult',1,'int',[1 1],'out',1,'in',1), struct('mult',2,'int',[2 2],'out',2,'in',2)};
pf_qdirect('S3_edge1d_heat2d',{struct('r',0.9*pi^2)},qdegs,[0 1],'s3');
pf_qdirect('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},qdegs,[0 1],'s17');
pf_qdirect('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},qdegs,[0 1],'s18');
disp('## QD RUN DONE');
