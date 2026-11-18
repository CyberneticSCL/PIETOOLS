% PF_RUN_P3_CAP2  Phase 3, the capped all-pattern row (Rlight_cap1: 2-D entry
% all, d1 1, d2 2, d3 [1 1;2 2;2 2], per-direction cap 1) as the light 2-D
% entry candidate: the rest of the mixed suite (S3 and S18 are in p3cap),
% light tier. The 2-D library rows follow in a second MATLAB call
% (PF_RUN_P2_LIB with STO = {'Rlight_cap1'}). SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
a26 = pi^2/4 - 0.5;
L = {'S1','S1_ode_heat1d_bc',struct('a',0.5),'stabP';
     'S8','S8_two_edges',struct('r3',0.9*pi^2),'stabP';
     'S4','S4_wentzell',struct('r',0.9*(pi^2+0.860333589^2)),'stabP';
     'S6','S6_indomain_1d_2d',struct('g',0.9*pi^2),'stabP';
     'S7','S7_ode_2d_indomain',struct('g',0.9),'stabP';
     'S9','S9_ode_2d_flux',struct('a',0.9),'stabP';
     'S26','S26_lib42_ode_heat_delay',struct('a',0.9*a26),'stabP';
     'S2','S2_ode_heat2d_edge',struct('a',0.5),'stabP';
     'S5','S5_chain_ode_1d_2d',struct('a',0.5),'stabP';
     'S23','S23_hinf_indomain',struct(),'hinfP';
     'S24','S24_hinf_ode_2d_edge',struct(),'hinfP';
     'S25','S25_hinf_wentzell',struct(),'hinfP'};
out = sprintf('%s/results/p3_cap2.mat',S);   Out = [];
for i = 1:size(L,1)
    try
        Sys = mixed_domain_systems(L{i,2},L{i,3});
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
    catch ME
        fprintf('@@ EX %s build ERROR %s\n',L{i,1},strrep(ME.message,newline,' '));  continue
    end
    r = pf_exec_row(L{i,1},PIE,L{i,4},'Rlight_cap1','light');
    if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end             %#ok<SAGROW>
    save(out,'Out');
end
disp('## P3CAP2 DONE');
