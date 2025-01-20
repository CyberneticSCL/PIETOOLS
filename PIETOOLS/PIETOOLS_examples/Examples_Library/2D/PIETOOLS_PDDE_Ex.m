function [PDE_t] = PIETOOLS_PDDE_Ex(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'lam=1;','dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based 
%               format.
%
% %---------------------------------------------------------------------% %
% Original 1D PDE: x1_{t}(t,s1) = x1_{s1s1}(t,s1) - 20 x1(t,s1) - 4 x4(t,s1,-1) - 0.1x5(t,s1,-1) + s1(s1-\pi)w(t);
% 2D Transport Equation resulting from delay 1: x4_{t}(t,s1,s2) = 4 x4_{s2}((t,s1,s2);
%  2D Transport Equation resulting from delay 2: x5_{t}(t,s1,s2) = 3.3333 x5_{s2}(t,s1,s2);
%  output z(t) = int_0^1 x1(t,s1) ds1;
%  BC1 0 = x1(t,0);
%  BC2 0 = x1(t,3.1416);
% BC3 0 = x4(t,s1,0) - x1(t,s1);
% BC4  0 = x5(t,s1,0) - x1(t,s1);
% H2 norm upper-bounds computed bwith custom settings and
% settings.settings_2d.ineq_opts.psatz=1: 0.8502(primal) and 0.8398(dual)
% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','H2_norm = 1;');
evalin('base','H2_norm_dual = 1;');

npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% Declare independent variables (time and space)
pvar s t
tau = [0.25;0.3];   
x1 = pde_var('state',1,s1,[0,pi]);
u = pde_var('input',1,[],[]);
y= pde_var('output',1,[],[]);
x2 = pde_var('state',1,[s1 s2],[0,pi;-1,0]);
x3 = pde_var('state',1,[s1 s2],[0,pi;-1,0]);
PDE1=diff(x1,'t')==diff(x1,s1,2)-20*x1-4*subs(x2,s2,-1)-0.1*subs(x3,s2,-1)+(s1*(s1-pi))*u;
PDE2=diff(x2,'t')==(1/tau(1))*diff(x2,s2);
PDE3=diff(x3,'t')==(1/tau(2))*diff(x3,s2);
BC1=subs(x1,s1,0)==0;
BC2=subs(x1,s1,pi)==0;
BC3=subs(x2,s2,0)==x1;
BC4=subs(x3,s2,0)==x1;
output=y==int(x1,s1,[0,pi]);
PDE_t=[PDE1;PDE2;PDE3;BC1;BC2;BC3;BC4;output];


if GUI~=0
    disp('No GUI representation available for this system.')
end

end