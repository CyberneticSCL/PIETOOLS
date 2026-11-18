function [PDE_t, st] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman(GUI,params) % MMP, 11/17/2026
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
% % Example 2), first-order (Riemann) form on the half interval:
% % PDE         z_{tt}(t,s) = 0.1*z_{ss} - 2*z_{t} - z(t-tau,s);
% %                                                       s in [0,pi]
% % With BCs    z(t,s=0) = 0;   z(t,s=pi) = 0;
% %
% % Restricted to the modes symmetric about s = pi/2 (the odd sine modes,
% % which contain the critical mode k = 1): s1 in [0,pi/2], z(t,0) = 0,
% % z_{s}(t,pi/2) = 0. Let e = z_{s}, v = z_{t}, and
% % y(t,s1,s2) = e(t-s2,s1), s2 in [0,tau], so that
% % z(t-tau,s1) = int_0^s1 y(t,r,tau) dr:
% % PDE         e_{t} = v_{s1};
% %             v_{t} = 0.1*e_{s1} - 2*v - int_0^s1 y(t,r,s2=tau) dr;
% %             y_{t} = -y_{s2};
% % With BCs    v(t,s1=0) = 0;   e(t,s1=pi/2) = 0;   y(t,s1,s2=0) = e(t,s1).
% %
% % Parameter tau can be set; default tau = 0.9*tau*.
% %
% % Threshold [derived from the modes]: stable iff tau < tau* = 2.881401,
% % the delay margin of mode k = 1 (see
% % PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman).
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026): infeasible at tau = 0.1, 0.5, 0.9 and 1.1 times
% % tau*, with the stock light settings (840,208 decision variables) and
% % at 0.5, 0.9, 1.1 times tau* with the proposed light row of the
% % settings proposal of 11/15/2026 (844,033); no certificate in this
% % form at this tier.
% %                                                                         % MMP, 11/17/2026
% % Application-specific settings (output st): lpisettings_sop('light')     % MMP, 11/17/2026
% % with the free operator Q of the direct form at degree 2 in two          % MMP, 11/17/2026
% % variables, st.lpi.var.deg{2} = 2, for                                   % MMP, 11/17/2026
% % PIETOOLS_stability_sop(PIE,st,'Q'). Measured (gap study 11/17/2026,     % MMP, 11/17/2026
% % gap_results.md G4, MOSEK): tau = 0.9*tau* certified at 281,497          % MMP, 11/17/2026
% % decision variables without weighted slack terms and at 302,957 with     % MMP, 11/17/2026
% % the faces at offset 1 and wcross 0 (the light slack); 1.1*tau*: no      % MMP, 11/17/2026
% % certificate. The storage form 'P' certifies this example only at        % MMP, 11/17/2026
% % 0.1*tau* (2-D storage at the heavy lengths, st.sop.dw_feas(2) = 0:      % MMP, 11/17/2026
% % 532,276); infeasible at 0.1, 0.5, 0.9 and 1.1 times tau* with the       % MMP, 11/17/2026
% % light storage.                                                          % MMP, 11/17/2026
% % Through the executive with this st (MOSEK, 11/17/2026): certified at    % MMP, 11/17/2026
% % the default 0.9*tau*, 302,957 variables, no raise.                      % MMP, 11/17/2026
% %                                                                         % MMP, 11/18/2026
% % Example-specific settings of 11/18/2026 (output st), form 'Q',          % MMP, 11/18/2026
% % library default tau = 0.9*2.881401: lpisettings_sop('light') with       % MMP, 11/18/2026
% % the slack lift degree at its minimum, st.lpi.ineq.dD = 0, and Q of      % MMP, 11/18/2026
% % degree 2 in two variables, st.sop.deg_Q{2} = 2. MOSEK, 11/18/2026:      % MMP, 11/18/2026
% % certified at 225,705 decision variables, 6,608 rows, no raise           % MMP, 11/18/2026
% % (OPT_STATUS 1.00, eta 2.4e-11); the 11/17 st (deg 2) and heavy:         % MMP, 11/18/2026
% % 302,957; light: infeasible (563,813 after one raise); light in form     % MMP, 11/18/2026
% % 'P': infeasible (638,505).                                              % MMP, 11/18/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Second output st, the settings that certify the
%                   example (maintainer: application-specific settings
%                   where the tiers do not); the light tier certifies it
%                   in neither form (header);
% function [PDE_t] = PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman(GUI,params) % MMP, 11/17/2026 (was)
% MMP, 11/18/2026: st: slack lift degree at its minimum and Q of degree
%                   2 set through deg_Q, which the direct form reads in
%                   place of lpi.var.deg: 25% fewer variables (header);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2 r1

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
L = pi/2;
e = pde_var(s1,[0,L],1);
v = pde_var(s1,[0,L],1);
y = pde_var([s1;s2],[0,L;0,tau],[0;1]);         % no derivative in s1, 1st in s2
PDE_t = [diff(e,'t')==diff(v,s1);
         diff(v,'t')==0.1*diff(e,s1)-2*v-int(subs(subs(y,s2,tau),s1,r1),r1,[0,s1]);
         diff(y,'t')==-diff(y,s2);
         subs(v,s1,0)==0;           subs(e,s1,L)==0;
         subs(y,s2,0)==e];


% Application-specific settings (header): the light tier with Q of the      % MMP, 11/18/2026 (was); MMP, 11/17/2026
% direct form at degree 2 in two variables; built only when requested       % MMP, 11/18/2026 (was); MMP, 11/17/2026
% Example-specific settings (header): light, slack lift degree at its       % MMP, 11/18/2026
% minimum (dD = 0), Q of degree 2 in two variables (st.sop.deg_Q)           % MMP, 11/18/2026
if nargout>1                                                                % MMP, 11/17/2026
    st = lpisettings_sop('light');                                          % MMP, 11/17/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
    st.sop.deg_Q{2} = 2;                                                    % MMP, 11/18/2026
%    st.lpi.var.deg{2} = 2;                                                 % MMP, 11/18/2026 (was); MMP, 11/17/2026
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
