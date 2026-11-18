% PF_RUN_M2  Solver pairing (plan rule: an UNKNOWN exit is paired with a
% SeDuMi solve before it enters a table): S17 and S18 as in M1, the N-D
% default storage, facesS and box at dw 1, solver SeDuMi. SCRATCH, MMP,
% 11/16/2026.
C = pf_cfg({'facesS','box'},2);
tau43 = 1.034726;
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,1,{'def'},'light','m2_s17','sedumi');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,1,{'def'},'light','m2_s18','sedumi');
disp('## M2 DONE');
