function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise(GUI,params)
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
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), mu0 = 1, tau = 1: certified at mu1 = 0.5 and
% % mu1 = 2 with the stock light settings (21,024 decision variables) and
% % with the proposed light row of the settings proposal of 11/15/2026
% % (22,480).
% %                                                                         % MMP, 11/17/2026
% % No second output st: the light tier certifies this example              % MMP, 11/17/2026
% % (record above). Its storage is the proposed light row since             % MMP, 11/17/2026
% % 11/17/2026; its slack of 11/17/2026 (wcross 0, 2-D faces at offset      % MMP, 11/17/2026
% % 1) was not rerun.                                                       % MMP, 11/17/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier certifies the example, so
%                   no second output st (application-specific settings);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');

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
