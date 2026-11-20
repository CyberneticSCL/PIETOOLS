function [PDE_t,st] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control(GUI,params)  % MMP, 11/19/2026
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
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Example reaction-diffusion equation with boundary control through ODE
% % ODE         x_{t}(t) = u(t);
% % PDE       v_{t}(t,s) = lam*v(t,s) + v_{ss}(t,s) + s*(2-s)*w(t)
% % With BCs    v(t,s=0) = 0
% %             v_{s}(t,s=1) = x(t)
% % Output      z(t) = [x(t); int_{0}^{1}v(t,s)ds]
% % Parameter lam can be set, defaults to 5.
% % Optimal state feedback can achieve L2 gain bound of 4.99
% % Adapted from Example 22 from Shivakumar, 2022 (see bottom of file).
% %---------------------------------------------------------------------% %
%
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
% function [PDE_t] = PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control(GUI,params)  % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s

%%% Executive Function:
evalin('base','Hinf_control = 1;')

lam = 5;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input format
clear stateNameGenerator
pde_var state x v  input w  output z1 z2  control u
v.vars = s;     v.dom = [0,1];
PDE_t = [diff(x,'t')==u;
         diff(v,'t')==lam*v+diff(v,s,2)+s*(2-s)*w;
         z1==x;
         z2==int(v,s,[0,1]);
         subs(v,s,0)==0;
         subs(diff(v,s),s,1)==x];


if GUI
    disp('No GUI representation available for this system.')
end

% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('light');                                              % MMP, 11/19/2026
st.lpi.pos.len_com{1} = 1;                                                  % MMP, 11/19/2026
st.lpi.pos.len_sos{1} = [3;2;1];                                            % MMP, 11/19/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/19/2026
st.lpi.var.deg{1} = [0 2 2];                                                % MMP, 11/19/2026

end
% @article{shivakumar2022h_,
%   title={$ H\_ $\{$$\backslash$infty$\}$ $-optimal control of coupled ODE-PDE systems using PIE framework and LPIs},
%   author={Shivakumar, Sachin and Das, Amritam and Weiland, Siep and Peet, Matthew},
%   journal={arXiv preprint arXiv:2208.13104},
%   year={2022}
% }