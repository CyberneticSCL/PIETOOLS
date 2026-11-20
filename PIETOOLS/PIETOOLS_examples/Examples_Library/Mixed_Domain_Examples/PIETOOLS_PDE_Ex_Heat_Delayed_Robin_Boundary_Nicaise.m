function [PDE_t,st] = PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'mu1=2;','tau=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % Heat equation with a delayed Robin boundary term (Nicaise, Valein,
% % Fridman 2009, system (1), constant delay, a = 1):
% % PDE         u_{t}(t,s1) = u_{s1s1}(t,s1);                s1 in [0,pi]
% % With BCs    u(t,s1=0) = 0;
% %             u_{s1}(t,pi) = -mu0*u(t,pi) - mu1*u(t-tau,pi);
% %
% % Let y(t,s2) = u(t-tau*s2,pi), s2 in [0,1] (system (6) of the paper):
% % PDE         u_{t} = u_{s1s1};   y_{t} = -(1/tau)*y_{s2};
% % With BCs    u(t,s1=0) = 0;
% %             u_{s1}(t,s1=pi) = -mu0*u(t,pi) - mu1*y(t,s2=1);
% %             y(t,s2=0) = u(t,s1=pi);
% %
% % Parameters mu0, mu1 and tau can be set; default mu0 = 1, mu1 = 0.5,
% % tau = 1 (values chosen here; the paper gives none).
% %
% % The paper gives sufficient conditions for exponential stability only
% % [abstract]; no threshold is known here. Spectral abscissa by the method
% % of lines (heat finite differences with a ghost node, upwind transport,
% % N heat nodes, M transport nodes) [finite differences]:
% %   mu1 = 0.5: -0.762312 (N = 100, M = 200), -0.762261 (N = 200, M = 400);
% %   mu1 = 2:   -0.206656 (N = 100, M = 200), -0.198679 (N = 200, M = 400).
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier certifies the example, so
%                   no second output st (application-specific settings);
% MMP, 11/19/2026: Output st, the example-specific settings, and the form
%                   of its tests named in the workspace (stability_form,
%                   Hinf_gain_form, H2_norm_form), from the search of
%                   11/19/2026 (maintainer: no default settings);
%                   the entry of 11/17/2026 (no second output st) no longer
%                   holds;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% MMP, 11/19/2026: the dated test records of the header removed: they
%                   recorded tests and files the user does not have, not
%                   the problem;
% function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise(GUI,params)  % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026

mu0 = 1;    mu1 = 0.5;  tau = 1;
npars = length(params);
if npars~=0
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end


% % % Construct the PDE.
%%% pde_var input
clear stateNameGenerator
u = pde_var(s1,[0,pi]);
y = pde_var(s2,[0,1],1);
PDE_t = [diff(u,'t')==diff(u,s1,2);
         diff(y,'t')==-(1/tau)*diff(y,s2);
         subs(u,s1,0)==0;
         subs(diff(u,s1),s1,pi)==-mu0*subs(u,s1,pi)-mu1*subs(y,s2,1);
         subs(y,s2,0)==subs(u,s1,pi)];


% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('stripped');                                           % MMP, 11/19/2026
st.lpi.ineq.dD = 0;                                                         % MMP, 11/19/2026
st.lpi.ineq.psatz{2} = 'none';                                              % MMP, 11/19/2026
st.lpi.ineq.cap_inc = 'auto';                                               % MMP, 11/19/2026
st.sop.deg_Q{2} = [2 1 1];                                                  % MMP, 11/19/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{nicaise2009stability,
%   title={Stability of the heat and of the wave equations with boundary time-varying delays},
%   author={Nicaise, Serge and Valein, Julie and Fridman, Emilia},
%   journal={Discrete and Continuous Dynamical Systems - Series S},
%   volume={2},
%   number={3},
%   pages={559--581},
%   year={2009},
%   doi={10.3934/dcdss.2009.2.559}
% }
