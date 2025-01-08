function [PDE_t] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control(GUI,params)
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
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
%
% %---------------------------------------------------------------------% %
%  PDE:          x_{t} = 3x+(s^2+0.2)x_{ss} +s^2w(t)
%  With BC:     x(s=0) = 0, x_{s}(s=1)=0
%  And outputs:   z(t) = [int(x(t,s),s,0,1)
%                                            u]
% % 
% % Parameter c can be set.
% %---------------------------------------------------------------------% %

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','H2_control = 1;');

% Construct the PDE.
pvar s t;
x = pde_var('state',1,s,[0,1]);
z = pde_var('output',2);
w = pde_var('input',1);
u = pde_var('control',1);
%y=pde_var('sense',1);
PDE_t = [diff(x,t)==3*x+(s^2+0.2)*diff(x,s,2)+s*u+s^2*w; 
         subs(x,s,0)==0; 
         subs(diff(x,s),s,1)==0;
         z==[int(x,s,[0,1]); u]];
%y==subs(x,s,1)
if GUI
    %%% Associated GUI save file
        error("No GUI file associated with this example");
end

end