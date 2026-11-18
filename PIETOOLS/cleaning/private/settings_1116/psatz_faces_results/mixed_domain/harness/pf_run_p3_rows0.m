% PF_RUN_P3_ROWS0  Phase 3, the separable light row (Rlight0: 2-D entry
% separable, d1 0, d2 2, d3 [1 1;2 2], no cap = today's light lengths with
% the mixed patterns) on the 2-D ladder (fin2 gain, rd2 at 0.9) and on the
% mixed suite, light tier: the light 2-D entry candidate beside the uniform
% row (Rdef, P2) and the capped all-pattern row (Rlight_cap1, P3). SCRATCH,
% MMP, 11/16/2026.
S = pf_setup(false);
out = sprintf('%s/results/p3_rows0.mat',S);   Out = [];
pvar s1 s2 t
x = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);  z = pde_var('out',1);
PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+0.5*2*pi^2*x+w; z==int(x,[s1;s2],[0,1;0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
evalc('PIEg = initialize(convert(PDE,{''pie'',''silent''}));');
clear stateNameGenerator
x9 = pde_var(1,[s1;s2],[0,1;0,1]);
PDE9 = [diff(x9,t)==diff(x9,s1,2)+diff(x9,s2,2)+0.9*2*pi^2*x9;
        subs(x9,s1,0)==0; subs(x9,s1,1)==0; subs(x9,s2,0)==0; subs(x9,s2,1)==0];
evalc('PIE9 = initialize(convert(PDE9,{''pie'',''silent''}));');
r = pf_exec_row('fin2',PIEg,'hinfP','Rlight0','light');
fprintf('@@ EXGAP fin2 storage=Rlight0 tier=light gap=%+.3e\n',r.gam/0.068606888-1);
Out = r;
r = pf_exec_row('rd2_09',PIE9,'stabP','Rlight0','light');
Out(end+1) = r;     save(out,'Out');
tau43 = 1.034726;   a26 = pi^2/4 - 0.5;
L = {'S1','S1_ode_heat1d_bc',struct('a',0.5),'stabP';
     'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2),'stabP';
     'S8','S8_two_edges',struct('r3',0.9*pi^2),'stabP';
     'S4','S4_wentzell',struct('r',0.9*(pi^2+0.860333589^2)),'stabP';
     'S6','S6_indomain_1d_2d',struct('g',0.9*pi^2),'stabP';
     'S7','S7_ode_2d_indomain',struct('g',0.9),'stabP';
     'S9','S9_ode_2d_flux',struct('a',0.9),'stabP';
     'S26','S26_lib42_ode_heat_delay',struct('a',0.9*a26),'stabP';
     'S2','S2_ode_heat2d_edge',struct('a',0.5),'stabP';
     'S5','S5_chain_ode_1d_2d',struct('a',0.5),'stabP';
     'S18_05','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1),'stabP';
     'S23','S23_hinf_indomain',struct(),'hinfP';
     'S24','S24_hinf_ode_2d_edge',struct(),'hinfP';
     'S25','S25_hinf_wentzell',struct(),'hinfP'};
for i = 1:size(L,1)
    try
        Sys = mixed_domain_systems(L{i,2},L{i,3});
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
    catch ME
        fprintf('@@ EX %s build ERROR %s\n',L{i,1},strrep(ME.message,newline,' '));  continue
    end
    r = pf_exec_row(L{i,1},PIE,L{i,4},'Rlight0','light');
    Out(end+1) = r;                                                         %#ok<SAGROW>
    save(out,'Out');
end
disp('## P3ROWS0 DONE');
