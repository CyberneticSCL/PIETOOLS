% PF_RUN_P3_CAP  Phase 3, the per-direction cap knob: the 2-D ladder gain
% fin2 (exact 0.068606888), the estimator of library PDE 48 (stock light
% 0.0061794), rd2 at 0.9 and the mixed rows S3 (0.9 r*) and S18 (0.5 tau*)
% under Rlight, Rlight_cap2 (control: no restriction at these lengths) and
% Rlight_cap1 (theta_k + s_k <= 1 per direction), with stock for reference;
% light. SCRATCH, MMP, 11/16/2026.
S = pf_setup(false);
STO = {'stock','Rlight','Rlight_cap2','Rlight_cap1'};
out = sprintf('%s/results/p3_cap.mat',S);   Out = [];
pvar s1 s2 t
x = pde_var(1,[s1;s2],[0,1;0,1]);  w = pde_var('in',1);  z = pde_var('out',1);
PDE = [diff(x,t)==diff(x,s1,2)+diff(x,s2,2)+0.5*2*pi^2*x+w; z==int(x,[s1;s2],[0,1;0,1]);
       subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
evalc('PIEg = initialize(convert(PDE,{''pie'',''silent''}));');
for sto = STO
    r = pf_exec_row('fin2',PIEg,'hinfP',sto{1},'light');
    fprintf('@@ EXGAP fin2 storage=%s gap=%+.3e\n',sto{1},r.gam/0.068606888-1);
    if isempty(Out),    Out = r;  else,   Out(end+1) = r;   end             %#ok<SAGROW>
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
for sto = STO
    r = pf_exec_row('PDE48',PIE48,'est',sto{1},'light');
    Out(end+1) = r;                                                         %#ok<SAGROW>
    save(out,'Out');
end
tau43 = 1.034726;
L = {'S3','S3_edge1d_heat2d',struct('r',0.9*pi^2); 'S18_05','S18_lib43_interior_delay',struct('tau',0.5*tau43,'sgn',-1)};
for i = 1:size(L,1)
    Sys = mixed_domain_systems(L{i,2},L{i,3});
    PIE = initialize(convert(Sys.PDE,{'pie','silent'}));
    for sto = STO
        r = pf_exec_row(L{i,1},PIE,'stabP',sto{1},'light');
        Out(end+1) = r;                                                     %#ok<SAGROW>
        save(out,'Out');
    end
end
disp('## P3CAP DONE');
