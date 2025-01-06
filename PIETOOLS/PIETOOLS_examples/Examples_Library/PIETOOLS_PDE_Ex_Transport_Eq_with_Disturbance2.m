function [PDE_t] = PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
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
%  PDE :        x_{t} = x_{s} + (s-s^2)w(t)
%  With BC     x(s=1) = 0
%  And output  z(t) = int(x(t,s),s,0,1)
% H2 norm computed by numerical integration of the definition: 0.1016
% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pde_struct PDE_t;
pvar s1 s2

%%% Executive Function:
evalin('base','H2_norm = 1;');
evalin('base','H2_norm_dual = 1;');

nu = 1;     r = 8;   ne = 1;
a = 0;  b=1;        c = 0;       d = 1;
Cw = (s1^2-b)*(s2^2-d);
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end

% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = pde_var('state',1,s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
% Declare the sytem equations
pde = [diff(x,t,1)==diff(x,s,1)+(s-s^2)*w;    % dynamics
                z==int(x,s,[0,1]);                     % output equation
                subs(x,s,1)==0];                            % boundary condition
PDE_t=initialize(pde);


if GUI~=0
    disp('No GUI representation available for this system.')
end

end