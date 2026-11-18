% PF_RUN_P3_HEAVY  Phase 3, the heavy 2-D entry: today's heavy lengths as a
% row (Rheavy_u: uniform, d3a 2; Rheavy_s: the same with the separable
% patterns) against the stock heavy storage on the 2-D ladder (fin2 gain,
% rd2 at 0.9; heavy tier), on the PDE 48 estimator (light tier, stock at
% heavy for reference) and, Rheavy_s, on S3 (0.9 r*) and S18 (0.5 tau*);
% then, if the clock allows, the 3-D heat at 0.9 kappa* under the 3-D
% per-direction cap 1 (Rlight_cap3). SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
out = sprintf('%s/results/p3_heavy.mat',S);   Out = [];
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
for sto = {'Rheavy_u','Rheavy_s'}
    r = pf_exec_row('fin2',PIEg,'hinfP',sto{1},'heavy');
    fprintf('@@ EXGAP fin2 storage=%s tier=heavy gap=%+.3e\n',sto{1},r.gam/0.068606888-1);
    if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end             %#ok<SAGROW>
    r = pf_exec_row('rd2_09',PIE9,'stabP',sto{1},'heavy');
    Out(end+1) = r;                                                         %#ok<SAGROW>
    save(out,'Out');
end
H = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad';
LB = [H '/psatz_mix/route4/lb'];    addpath(LB);
Cat = load([LB '/catalog.mat']);  C = Cat.C;
m = C(strcmp({C.lib},'PDE') & strcmp({C.id},'48'));    m = m(1);
params = {};
if m.idxarg,   PDE_t = feval(m.fn,48,0,params);   else,   PDE_t = feval(m.fn,0,params);   end
PDE48 = initialize_PIETOOLS_PDE(PDE_t,true);
evalc('PIE48 = initialize(convert(PDE48,{''pie'',''silent''}));');
L = {'stock','heavy'; 'Rheavy_u','light'; 'Rheavy_s','light'};
for i = 1:size(L,1)
    r = pf_exec_row('PDE48',PIE48,'est',L{i,1},L{i,2});
    Out(end+1) = r;                                                         %#ok<SAGROW>
    save(out,'Out');
end
tau43 = 1.034726;
L = {'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2); 'S18_05','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1)};
for i = 1:size(L,1)
    Sys = mixed_domain_systems(L{i,2},L{i,3});
    PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
    r = pf_exec_row(L{i,1},PIE,'stabP','Rheavy_s','light');
    Out(end+1) = r;                                                         %#ok<SAGROW>
    save(out,'Out');
end
% 3-D: the per-direction cap in three variables, if before 13:30
c = clock;
if c(4)*60 + c(5) < 13*60 + 30
    try
        pvar s3
        x3 = pde_var(1,[s1;s2;s3],[0,1;0,1;0,1]);
        PDE3 = [diff(x3,t)==diff(x3,s1,2)+diff(x3,s2,2)+diff(x3,s3,2)+0.9*1.5*pi^2*x3;
                subs(x3,s1,0)==0; subs(x3,s1,1)==0; subs(x3,s2,0)==0; subs(diff(x3,s2),s2,1)==0; subs(x3,s3,0)==0; subs(diff(x3,s3),s3,1)==0];
        evalc('PIE3 = initialize(convert(PDE3,{''pie'',''silent''}));');
        % 08:40 (D49): stock and Rlight were both infeasible at 0.9 kappa* under
        % the light 3-D settings (dw_feas 1, no raise, cap_delta 0: slack
        % bases bounded by Jmax(T)); rerun without that bound, all three rows
        for sto = {'stock','Rlight','Rlight_cap3'}
            r = pf_exec_row('heat3d_09_nocap',PIE3,'stabP',sto{1},'light',struct('cap_delta',[NaN NaN NaN]));
            Out(end+1) = r;                                                 %#ok<SAGROW>
            save(out,'Out');
        end
    catch ME
        fprintf('@@ EX heat3d_cap ERROR %s\n',strrep(ME.message,newline,' '));
    end
else
    disp('@@ EX heat3d_cap skipped (clock)');
end
disp('## P3HEAVY DONE');
