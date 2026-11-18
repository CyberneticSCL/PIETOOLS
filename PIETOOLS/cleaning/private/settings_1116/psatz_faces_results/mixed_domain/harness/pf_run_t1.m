% PF_RUN_T1  Phase 1, T1: the lift of the PDE norm gam I - T'*T by term set
% (PF_POINC, bracket [lam, 1.1 lam], 10 steps, light lengths) on S17, S3,
% S8, S18, S26. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
tau43 = 1.034726;   a26 = pi^2/4 - 0.5;
[optl,names] = pf_termsets();
sys = {'S17','S17_transport_edge_diffusion',struct('r2',0.5*pi^2); 'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2); ...
       'S8','S8_two_edges',struct('r3',0.9*pi^2); 'S18','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1); ...
       'S26','S26_lib42_ode_heat_delay',struct('a',0.9*a26)};
for i = 1:size(sys,1)
    try
        Sys = mixed_domain_systems(sys{i,2},sys{i,3});
        PIE = exec_pie_sop(initialize(convert(Sys.PDE,{'pie','silent'})));
        [sp,~] = copvar_space_list(PIE.T,'in');
        degs = zeros(1,numel(sp));
        for k = 1:numel(sp),    degs(k) = [0 6 4 3]*(numel(sp{k})==[0 1 2 3])';  end
        fprintf('@@ T1 %s terms=%s\n',sys{i,1},strjoin(names,','));
        pf_poinc(['t1_' sys{i,1}],sys{i,1},PIE,degs,optl,10,1.1);
    catch ME
        fprintf('@@ T1 ERROR %s: %s\n',sys{i,1},strrep(ME.message,newline,' '));
    end
end
disp('## T1 DONE');
