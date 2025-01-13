function [PDE_t] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_Control(GUI,params)
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
% % PDE                             x_{t} = 5*x + x_{ss} + w
% % With BCs                         x(s=0) = 0
% %                                         x_{s}(s=1) = xo
% % and regulated output  z =[ int(x(s,t),s,0,1);xo]
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
xo = pde_var('state',1,[],[]);
x = pde_var('state',1,s,[0,1]);
z = pde_var('output',2);
w = pde_var('input',1);
u = pde_var('control',1);
%y=pde_var('sense',1);
PDE_t = [diff(xo,t)==u;
                diff(x,t)==5*x+diff(x,s,2)+w; 
               subs(x,s,0)==0; 
              subs(diff(x,s),s,1)==xo;
              z==[int(x,s,[0,1]); xo]
              ];
%y==subs(x,s,1)
if GUI
    %%% Associated GUI save file
        error("No GUI file associated with this example");
end

end