% INH_RUN_3  Item 3: library PDE 7 (Gahlawat 2017) stabP light under stock
% and Rlight0 at lam 4.0, 4.5; bisection on [3,4.66] with stock at light and
% heavy; library PDE 23 (Shivakumar 2019 Ex. 1) hinfP stock light, heavy.
% SCRATCH, 11/16/2026.
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces');
S = pf_setup(false);
D = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/inhom';
addpath(D);
LB = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_mix/route4/lb';
addpath(LB);  Cat = load([LB '/catalog.mat']);  C = Cat.C;
m7 = C(strcmp({C.lib},'PDE') & strcmp({C.id},'7'));  m7 = m7(1);   disp(m7);
mk = @(lam) {sprintf('lam=%.10g;',lam),'c = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;'};
Out = [];   out = [D '/results/item3.mat'];
for lam = [4.0 4.5]
    for sto = {'stock','Rlight0'}
        PIE = lib7(m7,mk(lam));
        r = pf_exec_row(sprintf('G7_%.2f',lam),PIE,'stabP',sto{1},'light');  r.lam = lam;
        if isempty(Out),  Out = r;  else,  Out(end+1) = r;  end    %#ok<AGROW>
        save(out,'Out');
    end
end
for tier = {'light','heavy'}
    lo = 3;  hi = 4.66;  anyc = false;
    for it = 1:8
        lam = (lo+hi)/2;  PIE = lib7(m7,mk(lam));
        r = pf_exec_row(sprintf('G7_bis%d_%.5f',it,lam),PIE,'stabP','stock',tier{1});  r.lam = lam;
        Out(end+1) = r;  save(out,'Out');                                  %#ok<AGROW>
        if r.certified,  lo = lam;  anyc = true;  else,  hi = lam;  end
    end
    fprintf('@@ BIS G7 storage=stock %s lamcert=%.6f vs 4.66 reach=%.5f anycert=%d\n',tier{1},lo,lo/4.66,anyc);
end
m23 = C(strcmp({C.lib},'PDE') & strcmp({C.id},'23'));  m23 = m23(1);  disp(m23);
for tier = {'light','heavy'}
    PIE = lib7(m23,{});
    r = pf_exec_row('S19ex1_4.6',PIE,'hinfP','stock',tier{1});  r.lam = 4.6;
    Out(end+1) = r;  save(out,'Out');                                       %#ok<AGROW>
end
disp('## ITEM3 DONE');
function PIE = lib7(m,params)
clear stateNameGenerator
if m.idxarg,  PDE_t = feval(m.fn,str2double(m.id),0,params);  else,  PDE_t = feval(m.fn,0,params);  end
PDE = initialize_PIETOOLS_PDE(PDE_t,true);
evalc('PIE = initialize(convert(PDE,{''pie'',''silent''}));');
end
