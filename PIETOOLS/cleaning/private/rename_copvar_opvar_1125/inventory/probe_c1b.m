% C1 probe: is transl_mod compatible with a 1-D container? transl_mod_c is transl_mod with the
% input check widened to 'copvar', which is what the check becomes after the rename. Each
% container result is compared kernel by kernel with the struct path that PIESIM uses today.
root = 'C:/Users/mpeet/OneDrive/Documents/GitHub/PIETOOLS/PIETOOLS';
addpath(genpath(root));
addpath('C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/874c85cb-c18d-4876-bea3-39cdd053e21f/scratchpad/ren1125');
fprintf('transl_mod   : %s\ntransl_mod_c : %s\n', which('transl_mod'), which('transl_mod_c2'));
lastwarn('');

% Test 1: a container with every 1-D kernel nonzero, on [0,2]
ops = {};
try
    pvar s s_dum
    copvar T
    T.I = [0 2];  T.var1 = s;  T.var2 = s_dum;
    T.P = 2;  T.Q1 = s+1;  T.Q2 = s^2-1;
    T.R.R0 = 3+s;  T.R.R1 = s*s_dum+1;  T.R.R2 = s-2*s_dum;
    ops(end+1,:) = {'built, all kernels on [0,2]', T};
catch e, fprintf('T1 build: error [%s] %s\n', e.identifier, e.message); end

% Test 2: T and A of a converted 1-D PDE on [0,2] (heat equation, Dirichlet)
try
    pvar t s
    x = pde_var(1,s,[0,2]);
    PDE = initialize([diff(x,t) == diff(x,s,2); subs(x,s,0) == 0; subs(x,s,2) == 0]);
    PIE = convert(PDE,'pie');
    ops(end+1,:) = {'heat PIE.T', PIE.T};
    ops(end+1,:) = {'heat PIE.A', PIE.A};
catch e, fprintf('T2 build: error [%s] %s\n', e.identifier, e.message); end

for k = 1:size(ops,1)
    T = ops{k,2};
    fprintf('\n== %s: class %s, is_opvar1d %d\n', ops{k,1}, class(T), is_opvar1d(T));
    try
        St = transl_mod(kernels_of(T),[-1 1]);
    catch e, fprintf('  struct path: error [%s] %s\n', e.identifier, e.message); continue, end
    try
        Tc = transl_mod_c2(T,[-1 1]);
    catch e, fprintf('  container path: error [%s] %s\n', e.identifier, e.message); continue, end
    fprintf('  container path ran; result class %s, I = %s\n', class(Tc), mat2str(Tc.I));
    S2 = transl_mod_c2(kernels_of(T),[-1 1]);
    Kc = kernels_of(Tc);
    f = {'P','Q1','Q2','R0','R1','R2'};
    for i = 1:numel(f)
        a = getk(Kc,f{i});  b = getk(St,f{i});
        fprintf('  %-3s max|container - struct| = %.3g   (max|struct| = %.3g)\n', f{i}, pmax(a-b), pmax(b));
    end
end
[msg,id] = lastwarn;  fprintf('\nlast warning: [%s] %s\n', id, msg);

function K = kernels_of(T)
% The struct of kernels that pie2piesim hands to transl_mod, read through the 1-D emulation.
K = struct('P',T.P,'Q1',T.Q1,'Q2',T.Q2,'R',struct('R0',T.R.R0,'R1',T.R.R1,'R2',T.R.R2), ...
           'I',T.I,'var1',T.var1,'var2',T.var2);
end
function v = getk(K,f)
if f(1)=='R' && numel(f)==2,  v = K.R.(f);  else,  v = K.(f);  end
end
function m = pmax(X)
if isa(X,'polynomial'),  c = X.coefficient;  else,  c = double(X);  end
if isempty(c),  m = 0;  else,  m = full(max(abs(c(:))));  end
end
