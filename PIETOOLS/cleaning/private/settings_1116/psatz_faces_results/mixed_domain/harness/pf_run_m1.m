% PF_RUN_M1  S17 (transport state as Dirichlet edge data, r2* = pi^2) at
% r2 = 0.5 r2* and S18 (library PDE 43, sign -1, the delay a 2-D transport
% state, tau* = 1.034726) at tau = 0.5 tau*: the executive path (stability
% 'P', light, MOSEK) with the N-D default storage, the term sets none, box
% (Putinar, per-direction offset 1), facesS (today's 'auto'), faces0 and
% corners, dw 0 and 1 pinned. SCRATCH, MMP, 11/16/2026.
C = pf_cfg({'none','box','facesS','faces0','corners'},2);
tau43 = 1.034726;
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,[0 1],{'def'},'light','m1_s17');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,[0 1],{'def'},'light','m1_s18');
disp('## M1 DONE');
