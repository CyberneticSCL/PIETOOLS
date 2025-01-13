function [PDE_t] =PIETOOLS_PDE_Ex_Wave_Eq_Boundary_Control(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'dom=[0,1]'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the pde_var structure.
%
% %---------------------------------------------------------------------% %
%% 
% % ODE                             xo_{t}= u
% % PDE                             x_{t} = [0 1;0 0]x +  [0 0;1 0] x_{ss}+ [0;1]w
% % With BCs                         [1 0]x(s=0) = 0
% %                                         [1 0]x_{s}(s=1) = xo
% % and regulated output  z =[ xo; [1 0]int(x(s,t),s,0,1)]
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);
 clear stateNameGenerator
% Initialize variables
pvar s t

%%% Executive Function:
evalin('base','Hinf_control = 1;')

% Construct the PDE.
% xo = pde_var('state',1,[],[]);
% x1 = pde_var(s,[0,1]);
% x2 = pde_var(s,[0,1]);
% z = pde_var('output',2);
% w = pde_var('input',1);
% u = pde_var('control',1);
% %y=pde_var('sense',1);
% PDE_t = [diff(xo,t)==u;
%          diff(x1,t)==x2
%          diff(x2,t)==diff(x1,s,2)+w; 
%          subs(x1,s,0)==0; 
%          subs(diff(x1,s),s,1)==xo;
%          z==[int(x1,s,[0,1]); xo]];
% Declare state, input, and output variables
phi = pde_var('state',2,s,[0,1]);   x = pde_var('state',1,[],[]);
w = pde_var('input',1);             r = pde_var('output',1);
u = pde_var('control');   
% Declare system parameters
c=1;  
% declare dynamic equation
eq_dyn = [diff(x,t,1)==u
          diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;1]*w];
% declare output equation
eq_out= r ==int([1 0]*phi,s,[0,1]);
% declare the boundary conditions
bc1 = [0 1]*subs(phi,s,0)==0;   
bc2 = [1 0]*subs(phi,s,1)==x;
% create the PDE system
PDE_t = [eq_dyn;eq_out;bc1;bc2];
%y==subs(x,s,1)
if GUI
    %%% Associated GUI save file
        error("No GUI file associated with this example");
end

end