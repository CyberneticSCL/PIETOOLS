% PF_RUN_QD2  Soundness of the direct PIE2PDE certificates: S17 at 1.1 r2*
% (unstable: must fail) and 0.9 r2*, and at 0.5 r2* with eppos 0.1 (a
% margin-scaled artefact would not survive); S18 at 1.1 tau* and 0.9 tau*
% and at 0.5 tau* with eppos 0.1; Q at degree 1, dw 0. SCRATCH, MMP, 11/16/2026.
tau43 = 1.034726;
qd1 = {struct('mult',1,'int',[1 1],'out',1,'in',1)};
pf_qdirect('S17_transport_edge_diffusion',{struct('r2',1.1*pi^2),struct('r2',0.9*pi^2)},qd1,0,'s17_margin');
pf_qdirect('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},qd1,0,'s17_ep01',0.1);
pf_qdirect('S18_lib43_interior_delay',{struct('tau',1.1*tau43,'sgn',-1),struct('tau',0.9*tau43,'sgn',-1)},qd1,0,'s18_margin');
pf_qdirect('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},qd1,0,'s18_ep01',0.1);
disp('## QD2 RUN DONE');
