function inh_2d_sys(nm)
% INH_2D_SYS(NM) item 1: stabP at light on system NM at 0.9 lam* under
% 'stock' and 'Rlight0', then an 8-step bisection on [0.5,1] lam* under
% each storage. lam* = FD threshold at N = 60. SCRATCH, 11/16/2026.
D = 'C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/inhom';
L = load([D '/results/thresholds.mat']);    ls = L.T.(nm)(2);
out = sprintf('%s/results/item1_%s.mat',D,nm);  Out = [];
STO = {'stock','Rlight0'};
for k = 1:2
    PIE = inh_sys(nm,0.9*ls);
    r = pf_exec_row(sprintf('%s_0.9',nm),PIE,'stabP',STO{k},'light');  r.lam = 0.9*ls;
    if isempty(Out),  Out = r;  else,  Out(end+1) = r;  end    %#ok<AGROW>
    save(out,'Out');
end
for k = 1:2
    lo = 0.5;   hi = 1;     any_cert = false;
    for it = 1:8
        f = (lo+hi)/2;  PIE = inh_sys(nm,f*ls);
        r = pf_exec_row(sprintf('%s_bis%d_%.5f',nm,it,f),PIE,'stabP',STO{k},'light');  r.lam = f*ls;
        Out(end+1) = r;     save(out,'Out');                %#ok<AGROW>
        if r.certified,  lo = f;  any_cert = true;  else,  hi = f;  end
    end
    fprintf('@@ BIS %s storage=%s light lamstar=%.6f lamcert=%.6f reach=%.5f anycert=%d\n',nm,STO{k},ls,lo*ls,lo,any_cert);
end
fprintf('## ITEM1 %s DONE\n',nm);
end
