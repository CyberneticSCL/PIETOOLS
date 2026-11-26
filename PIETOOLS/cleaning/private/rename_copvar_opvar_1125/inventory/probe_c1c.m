% C1 check: the struct path (the one PIESIM uses) gives identical kernels through the original
% transl_mod and through the local-variable version transl_mod_c2.
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
addpath(genpath(root));
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/ren1125');
pvar s s_dum
copvar T
T.I = [0 2];  T.var1 = s;  T.var2 = s_dum;
T.P = 2;  T.Q1 = s+1;  T.Q2 = s^2-1;  T.R.R0 = 3+s;  T.R.R1 = s*s_dum+1;  T.R.R2 = s-2*s_dum;
pvar t
x = pde_var(1,s,[0,2]);
PIE = convert(initialize([diff(x,t) == diff(x,s,2); subs(x,s,0) == 0; subs(x,s,2) == 0]),'pie');
P2 = pie2piesim(PIE);                       % the structs PIESIM itself hands to rescalePIE
cases = {'built', kstruct(T); 'heat T (pie2piesim)', P2.T; 'heat A (pie2piesim)', P2.A};
for k = 1:size(cases,1)
    S0 = cases{k,2};
    A = transl_mod(S0,[-1 1]);   B = transl_mod_c2(S0,[-1 1]);
    fprintf('%-22s isequal(original, new) = %d\n', cases{k,1}, isequal(A,B));
end
function K = kstruct(T)
K = struct('P',T.P,'Q1',T.Q1,'Q2',T.Q2,'R',struct('R0',T.R.R0,'R1',T.R.R1,'R2',T.R.R2), ...
           'I',T.I,'var1',T.var1,'var2',T.var2);
end
