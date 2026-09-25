function [sol,M] = bl_b_stab1(phys,frac,setname,exec,n)
% 1-D stability. frac scales the destabilising term toward its analytic
% threshold, so 'how hard' is comparable across the three physics.
% n (optional) replicates the state n times -- decoupled, so the answer is
% unchanged and only the SDP grows. That is what makes the scaling arm a clean
% size sweep rather than a different problem at every rung.
if nargin<5 || isempty(n), n = 1; end
clear stateNameGenerator
pvar s t
switch phys
    case 'rd'                                  % lam* = pi^2, Dirichlet
        x = pde_var(n,s,[0,1]);
        PDE = [diff(x,t,1)==diff(x,s,2)+frac*pi^2*x; subs(x,s,0)==0; subs(x,s,1)==0];
    case 'tr'                                  % transport with damping
        x = pde_var(n,s,[0,1]);
        PDE = [diff(x,t,1)==-diff(x,s,1)-frac*x; subs(x,s,0)==0];
    case 'wave'                                % damped wave as a 2n-state system
        x = pde_var(2*n,s,[0,1]);
        PDE = [diff(x,t,1)==diff(x,s,2)-frac*x; subs(x,s,0)==0; subs(x,s,1)==0];
    otherwise, error('bl_b_stab1:phys','unknown physics %s',phys);
end
PIE = initialize(convert(PDE));
st  = bl_settings(setname);
[sol,P] = feval(['PIETOOLS_' exec],PIE,st);
M = struct('Pop',P,'Top',PIE.T,'Aop',PIE.A);
end


