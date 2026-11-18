function [PDE_t, st] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman(GUI,params) % MMP, 11/17/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'tau=1;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) settings structure that certifies the example    % MMP, 11/17/2026
%               (application-specific settings, settings report             % MMP, 11/17/2026
%               11/16/2026 Sec. 5.8).                                       % MMP, 11/17/2026
%
% %---------------------------------------------------------------------% %
% % Damped wave equation with an interior delay (Fridman, Orlov 2009,
% % Example 2, with the feedback u = -z(t-tau,s1)):
% % PDE         z_{tt}(t,s1) = 0.1*z_{s1s1} - 2*z_{t} - z(t-tau,s1);
% %                                                       s1 in [0,pi]
% % With BCs    z(t,s1=0) = 0;   z(t,s1=pi) = 0;
% %
% % Let v = z_{t}, and y(t,s1,s2) = z(t-s2,s1), s2 in [0,tau]:
% % PDE         z_{t} = v;
% %             v_{t} = 0.1*z_{s1s1} - 2*v - y(t,s1,s2=tau);
% %             y_{t} = -y_{s2};
% % With BCs    z = y = 0 at s1 = 0, pi;   y(t,s1,s2=0) = z(t,s1).
% % (v carries no derivative in s1.)
% %
% % Parameter tau can be set; default tau = 0.9*tau*.
% %
% % Threshold [derived from the modes]: mode k obeys
% % y'' + 2y' + 0.1*k^2*y + y(t-tau) = 0. A crossing of the imaginary axis
% % requires 0.1*k^2 < 1; the per-mode delay margins are 2.881401 (k = 1),
% % 3.497930 (k = 2), 8.758610 (k = 3), and modes k >= 4 are stable for
% % every delay:
% %   stable iff tau < tau* = 2.881401.
% % Rightmost root of mode 1: -0.1949 at 0.5*tau*, -0.0200 at 0.9*tau*,
% % +0.0155 at 1.1*tau*. The paper reports stability with decay rate 0.8
% % for tau <= 0.31 [read], a different question.
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), tau = 0.5*tau*: infeasible with the stock light
% % settings (691,888 decision variables) and with the proposed light row
% % of the settings proposal of 11/15/2026 (695,713).
% % PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman gives the
% % first-order (Riemann) form of the same system.
% %                                                                         % MMP, 11/17/2026
% % Application-specific settings (output st): lpisettings_sop('light')     % MMP, 11/17/2026
% % with the free operator Q of the direct form at degree 2 in two          % MMP, 11/17/2026
% % variables, st.lpi.var.deg{2} = 2, for                                   % MMP, 11/17/2026
% % PIETOOLS_stability_sop(PIE,st,'Q'). Measured (gap study 11/17/2026,     % MMP, 11/17/2026
% % gap_results.md G4, MOSEK): tau = 0.95*tau* certified at 235,353         % MMP, 11/17/2026
% % decision variables without weighted slack terms and at 248,497 with     % MMP, 11/17/2026
% % the faces at offset 1 and wcross 0 (the light slack); 1.1*tau*: no      % MMP, 11/17/2026
% % certificate. The storage form 'P' certifies this example at no          % MMP, 11/17/2026
% % setting measured: infeasible at 0.5*tau* and 0.95*tau* with the         % MMP, 11/17/2026
% % light storage, UNKNOWN at 0.95*tau* with the 2-D storage at the         % MMP, 11/17/2026
% % heavy or veryheavy lengths.                                             % MMP, 11/17/2026
% % Through the executive with this st (MOSEK, 11/17/2026): certified at    % MMP, 11/17/2026
% % 0.95*tau* and at the default 0.9*tau*, 248,497 variables, no raise.     % MMP, 11/17/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Second output st, the settings that certify the
%                   example (maintainer: application-specific settings
%                   where the tiers do not); the light tier certifies it
%                   in neither form (header);
% function [PDE_t] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman(GUI,params) % MMP, 11/17/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');

tau = 0.9*2.881401;
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
z = pde_var(s1,[0,pi]);
v = pde_var(s1,[0,pi],0);                       % v = z_{t}: no derivative in s1
y = pde_var([s1;s2],[0,pi;0,tau],[2;1]);        % 2nd order in s1, 1st in s2
PDE_t = [diff(z,'t')==v;
         diff(v,'t')==0.1*diff(z,s1,2)-2*v-subs(y,s2,tau);
         diff(y,'t')==-diff(y,s2);
         subs([z;y],s1,0)==0;       subs([z;y],s1,pi)==0;
         subs(y,s2,0)==z];


% Application-specific settings (header): the light tier with Q of the      % MMP, 11/17/2026
% direct form at degree 2 in two variables; built only when requested       % MMP, 11/17/2026
if nargout>1                                                                % MMP, 11/17/2026
    st = lpisettings_sop('light');                                          % MMP, 11/17/2026
    st.lpi.var.deg{2} = 2;                                                  % MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{fridman2009exponential,
%   title={Exponential stability of linear distributed parameter systems with time-varying delays},
%   author={Fridman, Emilia and Orlov, Yury},
%   journal={Automatica},
%   volume={45},
%   number={1},
%   pages={194--201},
%   year={2009},
%   doi={10.1016/j.automatica.2008.06.006}
% }
