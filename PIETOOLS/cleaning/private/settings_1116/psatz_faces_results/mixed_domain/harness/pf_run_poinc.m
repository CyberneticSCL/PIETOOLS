% PF_RUN_POINC  gam I - T'*T on the PIEs of S17, S3, S18, S17's x2 alone and
% the 1-D transport alone (PF_POINC): Rayleigh lam over Legendre products
% (degree 6 on a 1-D space, 4 per variable on a 2-D space), then the
% certified reach gam_min/lam under no terms, the reader's default terms,
% and raised degrees (8 bisection steps). SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
pvar s1 s2 t
tau43 = 1.034726;
optl = {struct('psatz','none'), struct(), struct('dD',2,'dw',2)};
sys = {'S17','S17_transport_edge_diffusion',struct('r2',0.5*pi^2); 'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2); 'S18','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1)};
for i = 1:size(sys,1)
    Sys = mixed_domain_systems(sys{i,2},sys{i,3});
    PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
    pf_poinc(['poinc_' sys{i,1}],sys{i,1},PIE,[6 4],optl,8);
end
clear stateNameGenerator
x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);
PDE = [diff(x2,t)==diff(x2,s2,2); subs(x2,s2,0)==0; subs(x2,s2,1)==0];
PIE = exec_pie_sop(initialize(convert(PDE,{'pie','silent'})));
pf_poinc('poinc_S17x2','S17x2',PIE,[4],optl,8);
clear stateNameGenerator
x1 = pde_var(1,s1,[0,1]);
PDE = [diff(x1,t)==-diff(x1,s1); subs(x1,s1,0)==0];
PIE = exec_pie_sop(initialize(convert(PDE,{'pie','silent'})));
pf_poinc('poinc_D1','D1_transport_1d',PIE,[6],optl,8);
disp('## POINC RUN DONE');
