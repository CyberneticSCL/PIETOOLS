% PF_RUN_M5  The decay-margin program (PF_DMAX) on S17 and S18 at the M1
% parameters: storages STO, term sets facesS and box, dw 1 pinned; delta* =
% 1 (the cap) is a strict certificate, 0 a program without interior.
% SCRATCH, MMP, 11/16/2026.
STO = {'def'};
C = pf_cfg({'facesS','box'},2);
tau43 = 1.034726;
pf_dmax('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,1,STO,'light','m5_s17');
pf_dmax('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,1,STO,'light','m5_s18');
disp('## M5 DONE');
