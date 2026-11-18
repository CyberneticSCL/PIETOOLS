% PF_RUN_P2_MIXED2  Phase 2, the rest of the mixed suite after PF_RUN_P2_MIXED
% was stopped at S9 (07:15): the separable row at d2 = 3 (Rlight2) costs 10x
% the stock storage (S7: 5.39M variables, 1940 s) and certifies nothing the
% lighter rows do not, except S18; so S26, S2, S5, S23, S24, S25 run under
% stock, Rdef, Rlight, and the two S18 rows under all four. Same tasks and
% margins as PF_RUN_P2_MIXED. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
tau43 = 1.034726;   a26 = pi^2/4 - 0.5;
STO3 = {'stock','Rdef','Rlight'};   STO4 = [STO3,{'Rlight2'}];
L = {'S26','S26_lib42_ode_heat_delay',struct('a',0.9*a26),'stabP',{STO3};
     'S2','S2_ode_heat2d_edge',struct('a',0.5),'stabP',{STO3};
     'S5','S5_chain_ode_1d_2d',struct('a',0.5),'stabP',{STO3};
     'S18_05','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1),'stabP',{STO4};
     'S18_09','S18_lib43_interior_delay',struct('tau',0.9*tau43,'sgn',-1),'stabP',{STO4};
     'S23','S23_hinf_indomain',struct(),'hinfP',{STO3};
     'S24','S24_hinf_ode_2d_edge',struct(),'hinfP',{STO3};
     'S25','S25_hinf_wentzell',struct(),'hinfP',{STO3}};
out = sprintf('%s/results/p2_mixed2.mat',S);   Out = [];
for i = 1:size(L,1)
    try
        Sys = mixed_domain_systems(L{i,2},L{i,3});
        PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
    catch ME
        fprintf('@@ EX %s build ERROR %s\n',L{i,1},strrep(ME.message,newline,' '));  continue
    end
    for sto = L{i,5}{1}
        r = pf_exec_row(L{i,1},PIE,L{i,4},sto{1},'light');
        if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end         %#ok<SAGROW>
        save(out,'Out');
    end
end
disp('## P2MIXED2 DONE');
