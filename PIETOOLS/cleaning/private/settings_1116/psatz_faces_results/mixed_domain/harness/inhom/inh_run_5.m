% INH_RUN_5  Follow-up: the mixed cascade M at 0.9 of both thresholds,
% stabP under stock at heavy (stock light returned UNKNOWN after one raise),
% then an 8-step bisection of the common fraction f (r1 = f r1*, r2 = f r2*)
% on [0.5,1] under stock and Rlight0 at light.
% SCRATCH, 11/16/2026.
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/psatz_faces');
S = pf_setup(false);
D = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/inhom';
addpath(D);
L = load([D '/results/thresholds.mat']);
r1 = 0.9*L.T.M1(2);   r2 = 0.9*L.T.M2(2);
PIE = inh_sys('M',[r1 r2]);
Out = pf_exec_row('M_0.9',PIE,'stabP','stock','heavy');
save([D '/results/item5.mat'],'Out');

for sto = {'stock','Rlight0'}
    lo = 0.5;  hi = 1;  anyc = false;
    for it = 1:8
        f = (lo+hi)/2;  PIE = inh_sys('M',f*[L.T.M1(2) L.T.M2(2)]);
        r = pf_exec_row(sprintf('M_bis%d_%.5f',it,f),PIE,'stabP',sto{1},'light');
        Out(end+1) = r;  save([D '/results/item5.mat'],'Out');             %#ok<SAGROW>
        if r.certified,  lo = f;  anyc = true;  else,  hi = f;  end
    end
    fprintf('@@ BIS M storage=%s light reach=%.5f anycert=%d\n',sto{1},lo,anyc);
end
disp('## ITEM5 DONE');
