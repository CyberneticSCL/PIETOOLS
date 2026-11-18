function [PDE_t, st] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params) % MMP, 11/17/2026
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
% - st:         (optional) settings structure that certifies the example    % MMP, 11/17/2026
%               (application-specific settings, settings report             % MMP, 11/17/2026
%               11/16/2026 Sec. 5.8).                                       % MMP, 11/17/2026
% - PDE_b:      PDE structure defining the example system in the batch
%               format.
%
% %---------------------------------------------------------------------% %
% % Diffusive Equation with Delay in Dynamics (Caliskan 2009):
% % PDE         x_{t}  = c*x_{s1s1}(t,s1) + a0*x(t,s1) - a1*x(t-tau,s1);    s1 in [0,pi]   
% % With BCs    x(t,s1=0) = 0;   
% %             x(t,s1=pi) = 0;
% %
% % Let x1=x, and x2(t,s1,s2)=x1(t-s2,s1), introducing transport equation:
% %         x1(t) = x2(t,0),   x2_{t} = -x2_{s2}        s2 in [0,tau]
% % Then
% %
% % PDE:    x1_{t}         = c*x1_{s1s1}(t,s1) + a0*x1(t,s1) - a1*x2(t,s1,s2=tau);
% %         x2_{t}         = -x2_{s2};
% % BCs:    x1(t,s1=0)     = 0;
% %         x1(t,s1=pi)    = 0;
% %         x2(t,s1,s2=0)  = x1(t,s1);
% %         x2(t,s1=0,s2)  = 0;
% %         x2(t,s1=pi,s2) = 0;
% %
% % Parameters tau, c, a0 and a1 can be set.
% % For c = 1, a0 = 1.9, and a1 = 1, stable for tau<1.0347.
% %                                                                         % MMP, 11/17/2026
% % Light tier (lpisettings_sop('light')), tau = 0.9*tau*: no certificate   % MMP, 11/17/2026
% % in the storage form 'P' (infeasible, 508,228 decision variables after   % MMP, 11/17/2026
% % one raise) or the direct form 'Q' (UNKNOWN, 227,510; gap_results.md     % MMP, 11/17/2026
% % G2, perspace_results.md, before the tier change of 11/17/2026).         % MMP, 11/17/2026
% % Application-specific settings (output st): lpisettings_sop('light')     % MMP, 11/17/2026
% % with the slack weight of the space of x1 (L2[s1]) raised by one,        % MMP, 11/17/2026
% % st.lpi.ineq.dw_space = [1 0], for PIETOOLS_stability_sop(PIE,st,'Q').   % MMP, 11/17/2026
% % Certified through the executive (MOSEK, 11/17/2026): tau = 0.9*tau*,    % MMP, 11/17/2026
% % 65,731 variables, no raise; default tau = 1 (0.966*tau*), 397,870       % MMP, 11/17/2026
% % after one raise. Alternative in 'Q': st.lpi.var.deg{2} = 2 (215,188     % MMP, 11/17/2026
% % at 0.9*tau*, gap study). In 'P': the 2-D storage at the veryheavy       % MMP, 11/17/2026
% % lengths, st.lpi.pos.len_com{2} = [3 3], st.lpi.pos.len_sos{2} =         % MMP, 11/17/2026
% % [2 2;3 3;0 0], with st.sop.dw_feas(2) = 0 (953,001 at 0.9*tau*, gap     % MMP, 11/17/2026
% % study, before the tier change of 11/17/2026).                           % MMP, 11/17/2026
% %                                                                         % MMP, 11/18/2026
% % Example-specific settings of 11/18/2026 (output st, replaces the        % MMP, 11/18/2026
% % dw_space setting above), form 'Q', library default tau = 1:             % MMP, 11/18/2026
% % lpisettings_sop('light') with the slack lift degree at its minimum,     % MMP, 11/18/2026
% % st.lpi.ineq.dD = 0, and Q of degree 2 in two variables,                 % MMP, 11/18/2026
% % st.sop.deg_Q{2} = 2. MOSEK, 11/18/2026: certified at 196,513            % MMP, 11/18/2026
% % decision variables, 6,376 rows, no raise (OPT_STATUS 0.998, eta         % MMP, 11/18/2026
% % 6.3e-10); the 11/17 st (dw_space = [1 0]): 397,870 after one raise      % MMP, 11/18/2026
% % (OPT_STATUS 1.16); heavy: 215,188, no raise; light: 374,827 after       % MMP, 11/18/2026
% % one raise. Light in form 'P': infeasible (428,185 after one raise).     % MMP, 11/18/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Second output st, the settings that certify the
%                   example (maintainer: application-specific settings
%                   where the tiers do not); the light tier certifies it
%                   in neither form (header);
% function [PDE_t] = PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay(GUI,params)   % MMP, 11/17/2026 (was)
% MMP, 11/18/2026: st: the slack lift degree at its minimum and Q of
%                   degree 2 in two variables (deg_Q), in place of the
%                   dw_space raise of 11/17: half the variables (header);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
%evalin('base','stability_dual = 1;');

c = 1;      a0 = 1.9;       a1 = 1;
ne = 1;     tau = 1; 
a = 0;      b = pi;
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
x1 = pde_var(s1,[a,b]);
x2 = pde_var([s1;s2],[a,b;0,tau],[2;1]);    % x2 must be 2nd order differentiable wrt s1, and 1st order wrt s2
PDE_t = [diff(x1,'t')==c*diff(x1,s1,2)+a0*x1-a1*subs(x2,s2,tau);
         diff(x2,'t')==-diff(x2,s2);
         subs([x1;x2],s1,a)==0;     subs([x1;x2],s1,b)==0;
         subs(x2,s2,0)==x1];


% %%% Term-based input format
% % Initialize the state variables and input.
% pde_struct PDE_t;
% PDE_t.x{1}.vars = s1;           PDE_t.x{1}.dom = [a,b];
% PDE_t.x{2}.vars = [s1;s2];      PDE_t.x{2}.dom = [a,b;0,tau];
% PDE_t.x{2}.diff = [2,1];    % x2 must be 2nd order differentiable wrt s1, and 1st order wrt s2
% 
% % PDE: x1_{t}(t,s1) = c*x1_{s1s1}(t,s1)
% PDE_t.x{1}.term{1}.x = 1;
% PDE_t.x{1}.term{1}.D = 2;
% PDE_t.x{1}.term{1}.C = c;
% 
% % PDE: x1_{t}(t,s1) = ... + a0*x1(t,s1)
% PDE_t.x{1}.term{2}.x = 1;
% PDE_t.x{1}.term{2}.C = a0;
% 
% % PDE: x1_{t}(t,s1) = ... + a1*x2(t,s1,tau)
% PDE_t.x{1}.term{3}.x = 2;
% PDE_t.x{1}.term{3}.loc = [s1,tau];
% PDE_t.x{1}.term{3}.C = -a1;
% 
% % PDE: x2_{t}(t,s1,s2) = -x2_{s2}(t,s1,s2)
% PDE_t.x{2}.term{1}.x = 2;
% PDE_t.x{2}.term{1}.D = [0,1];
% PDE_t.x{2}.term{1}.C = -1;
% 
% 
% % BC1: x1(s1=a) = 0;                BC2: x1(s1=b)=0;
% PDE_t.BC{1}.term{1}.x = 1;          PDE_t.BC{2}.term{1}.x = 1;
% PDE_t.BC{1}.term{1}.loc = a;        PDE_t.BC{2}.term{1}.loc = b;
% 
% % BC3: (Sewing condition) x2(s1,s2=0) = x1(s1)
% PDE_t.BC{3}.term{1}.x = 2;          PDE_t.BC{3}.term{2}.x = 1;
% PDE_t.BC{3}.term{1}.loc = [s1,0];   PDE_t.BC{3}.term{2}.C = -1;
% 
% % BC4: x2(s1=a,s2) = 0;             BC5: x2(s1=b,s2) = 0;
% PDE_t.BC{4}.term{1}.x = 2;          PDE_t.BC{5}.term{1}.x = 2;
% PDE_t.BC{4}.term{1}.loc = [a,s2];   PDE_t.BC{5}.term{1}.loc = [b,s2];


% Application-specific settings (header): the light tier, the slack         % MMP, 11/18/2026 (was); MMP, 11/17/2026
% weight of the space of x1 (L2[s1], the first space) raised by one         % MMP, 11/18/2026 (was); MMP, 11/17/2026
% Example-specific settings (header): light, slack lift degree at its       % MMP, 11/18/2026
% minimum (dD = 0), Q of degree 2 in two variables (st.sop.deg_Q)           % MMP, 11/18/2026
if nargout>1                                                                % MMP, 11/17/2026
    st = lpisettings_sop('light');                                          % MMP, 11/17/2026
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
    st.sop.deg_Q{2} = 2;                                                    % MMP, 11/18/2026
%    st.lpi.ineq.dw_space = [1 0];                                          % MMP, 11/18/2026 (was); MMP, 11/17/2026
end                                                                         % MMP, 11/17/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @article{CALISKAN2009,
% title = {Stability Analysis of the Heat Equation with Time-Delayed Feedback},
% journal = {IFAC Proceedings Volumes},
% volume = {42},
% number = {6},
% pages = {220-224},
% year = {2009},
% note = {6th IFAC Symposium on Robust Control Design},
% issn = {1474-6670},
% doi = {https://doi.org/10.3182/20090616-3-IL-2002.00038},
% url = {https://www.sciencedirect.com/science/article/pii/S1474667015404057},
% author = {Sina Yamaç çalişkan and Hitay özbay},
% }