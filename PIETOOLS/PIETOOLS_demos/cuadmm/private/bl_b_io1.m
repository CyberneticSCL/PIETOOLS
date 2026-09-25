function [sol,M] = bl_b_io1(exec,setname,n)
% 1-D plant WITH a disturbance and a regulated output -- the minimum structure
% the Hinf gain and H2 norm executives require. Without w the gain executives
% have no input channel and either error or return a degenerate program.
if nargin<3 || isempty(n), n = 1; end
clear stateNameGenerator
pvar s t
x = pde_var(n,s,[0,1]);
w = pde_var('in',n,s,[0,1]);
z = pde_var('out',n);
PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+w;
       z==int(x,s,[0,1]);
       subs(x,s,0)==0;  subs(x,s,1)==0];
PIE = initialize(convert(PDE));
st  = bl_settings(setname);
[sol,P,gam] = feval(['PIETOOLS_' exec],PIE,st);
M = struct('Pop',P,'Top',PIE.T,'Aop',PIE.A,'gam',gam);
end


