function [CC, ZdC] = minus_decparam(CA,CB,ZdA,ZdB)
minusCB = struct('A',-CB.A,'B',-CB.B);
[CC,ZdC] = plus_dpvar(CA,minusCB,ZdA,ZdB);
end