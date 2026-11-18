% PF_RUN_P2_MIXED  Phase 2, the mixed suite under the row storages: stability
% 'P' at 0.9 of the margin on S3, S8, S4, S6, S7, S9, S26 (S2, S5 at
% a = 0.5), S18 at 0.5 and 0.9 tau*; gains 'P' on S23, S24, S25 (exact
% 0.0093385782, 0.5, 0.0692045143); storages STO (default stock, Rdef,
% Rlight, Rlight2); light; the executive's own raise loop and terms.
% SCRATCH, MMP, 11/16/2026.
if ~exist('STO','var'),     STO = {'stock','Rdef','Rlight','Rlight2'};   end
S = pf_setup(false);
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
     'S18_09','S18_lib43_interior_delay',struct('tau',0.9*tau43,'sgn',-1),'stabP';
     'S23','S23_hinf_indomain',struct(),'hinfP';
     'S24','S24_hinf_ode_2d_edge',struct(),'hinfP';
     'S25','S25_hinf_wentzell',struct(),'hinfP'};
out = sprintf('%s/results/p2_mixed.mat',S);   Out = [];
for i = 1:size(L,1)
    try
        Sys = mixed_domain_systems(L{i,2},L{i,3});
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
    catch ME
        fprintf('@@ EX %s build ERROR %s\n',L{i,1},strrep(ME.message,newline,' '));  continue
    end
    for sto = STO
        r = pf_exec_row(L{i,1},PIE,L{i,4},sto{1},'light');
        if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end         %#ok<SAGROW>
        save(out,'Out');
    end
end
disp('## P2MIXED DONE');
