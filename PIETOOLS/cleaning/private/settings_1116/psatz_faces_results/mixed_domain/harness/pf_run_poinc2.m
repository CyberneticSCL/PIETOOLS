% PF_RUN_POINC2  The gap at the sharp constant: gam I - T'*T with the
% bracket [lam, 1.1 lam] and 12 bisection steps (resolution 2.4e-5 lam),
% the reader's default terms at light and at raised degrees, on S17, S3,
% S18 and the 1-D transport alone. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
pvar s1 s2 t
tau43 = 1.034726;
optl = {struct(), struct('dD',2,'dw',2)};
sys = {'S17','S17_transport_edge_diffusion',struct('r2',0.5*pi^2); 'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2); 'S18','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1)};
for i = 1:size(sys,1)
    Sys = mixed_domain_systems(sys{i,2},sys{i,3});
    PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
    pf_poinc(['poinc2_' sys{i,1}],sys{i,1},PIE,[6 4],optl,12,1.1);
end
clear stateNameGenerator
x1 = pde_var(1,s1,[0,1]);
PDE = [diff(x1,t)==-diff(x1,s1); subs(x1,s1,0)==0];
PIE = exec_pie_sop(initialize(convert(PDE,{'pie','silent'})));
pf_poinc('poinc2_D1','D1_transport_1d',PIE,[8],optl,12,1.1);
disp('## POINC2 RUN DONE');
