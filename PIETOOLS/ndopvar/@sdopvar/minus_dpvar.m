function [CC, ZdC] = minus_dpvar(CA,CB,ZdA,ZdB)
minusCB = struct('A',-CB.A,'Bt',-CB.Bt);
[CC,ZdC] = plus_dpvar(CA,minusCB,ZdA,ZdB);
end