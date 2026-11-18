% PF_RUN_M3  The storage question on S17 and S18 (same parameters as M1):
% every basis operator of every space at the light degrees, with the
% kernels split (lower and upper integrals) or separable (one full
% integral) per registry direction: 'sep11' both separable (the control for
% the extra mixed operators), 'sep01' s1 split, 'sep10' s2 split, 'sep00'
% both split. Term set facesS (today's 'auto'), dw 0 and 1 pinned, MOSEK.
% SCRATCH, MMP, 11/16/2026.
C = pf_cfg({'facesS'},2);
tau43 = 1.034726;
pf_exec_m('S17_transport_edge_diffusion',{struct('r2',0.5*pi^2)},C,[0 1],{'sep11','sep01','sep10','sep00'},'light','m3_s17');
pf_exec_m('S18_lib43_interior_delay',{struct('tau',0.5*tau43,'sgn',-1)},C,[0 1],{'sep11','sep10','sep01','sep00'},'light','m3_s18');
disp('## M3 DONE');
