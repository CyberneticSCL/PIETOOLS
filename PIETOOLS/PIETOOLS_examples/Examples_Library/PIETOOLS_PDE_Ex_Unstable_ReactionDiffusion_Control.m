function [PDE_t,st] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control(GUI,params)  % MMP, 11/19/2026
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
%  PDE:          x_{t} = 3x+(s^2+0.2)x_{ss} +s^2w(t)
%  With BC:     x(s=0) = 0, x_{s}(s=1)=0
%  And outputs:   z(t) = [int(x(t,s),s,0,1)
%                                            u]
% %---------------------------------------------------------------------% %
%
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t] = PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control(GUI,params)  % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s t

%%% Executive Function:
evalin('base','H2_control = 1;');

% Construct the PDE.
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

% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('heavy');                                              % MMP, 11/19/2026
st.lpi.pos.len_com{1} = 1;                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [3;2;1];                                            % MMP, 11/19/2026
st.lpi.var.deg{1} = [1 2 2];                                                % MMP, 11/19/2026
st.sop.dw_P = 1;                                                            % MMP, 11/19/2026
st.sop.deg_Q{1} = [6 2 2];                                                  % MMP, 11/19/2026

end