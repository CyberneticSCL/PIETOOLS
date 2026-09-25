function [sol,M] = bl_b_syn1(exec,setname)
% 1-D plant with BOTH a control input and a sensed output, so one plant serves
% the estimator, the controller and well-posedness. An estimator case on a plant
% with no sensed output, or a controller case with no control input, is not a
% harder case -- it is a different, degenerate program.
%
% Output count from nargout(), and M assembled field by field: struct() with a
% cell value expands into a struct ARRAY rather than storing the cell.
clear stateNameGenerator
pvar s t
x = pde_var(1,s,[0,1]);
w = pde_var('in',1,s,[0,1]);
u = pde_var('control',1);
y = pde_var('sense',1);
z = pde_var('out',2);
PDE = [diff(x,t,1)==diff(x,s,2)+0.5*pi^2*x+w+u;
       z==[int(x,s,[0,1]); u];
       y==int(x,s,[0,1]);
       subs(x,s,0)==0;  subs(x,s,1)==0];
PIE = initialize(convert(PDE));
st  = bl_settings(setname);
fn  = ['PIETOOLS_' exec];
nout = abs(nargout(fn));
out = cell(1,nout);
[out{1:nout}] = feval(fn,PIE,st);
sol = out{1};

M = struct();
M.Top = PIE.T;
M.Aop = PIE.A;
M.aux = {out(2:end)};
for k = 2:nout
    if isnumeric(out{k}) && isscalar(out{k}), M.gam = double(out{k}); break; end
end
end
