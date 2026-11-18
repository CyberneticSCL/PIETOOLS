function [PDE_t] = PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'r=2.9;','tau=0.5;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
%
% %---------------------------------------------------------------------% %
% % 2D heat equation with an interior delay (2D extension of Fridman,
% % Orlov 2009, Example 1):
% % PDE         x_{t}(t,s1,s2) = x_{s1s1} + x_{s2s2} + r*x(t,s1,s2)
% %                                - x(t-tau,s1,s2);    (s1,s2) in [0,pi]^2
% % With BCs    x = 0 at s1 = 0, s1 = pi, s2 = 0, s2 = pi.
% %
% % Let x1 = x, and x2(t,s1,s2,s3) = x1(t-s3,s1,s2), s3 in [0,tau]:
% % PDE         x1_{t} = x1_{s1s1} + x1_{s2s2} + r*x1 - x2(t,s1,s2,s3=tau);
% %             x2_{t} = -x2_{s3};
% % With BCs    x1 = x2 = 0 at s1 = 0, pi and at s2 = 0, pi;
% %             x2(t,s1,s2,s3=0) = x1(t,s1,s2).
% % The 2D state x1 is the datum of the 3D state x2 on the face s3 = 0.
% %
% % Parameters r and tau can be set; default r = 2.9, tau = 0.9*tau*.
% %
% % Threshold [derived from the modes, exact]: the sine mode (m,n) obeys
% % y' = a*y - y(t-tau), a = r-m^2-n^2, with delay margin
% % arccos(a)/sqrt(1-a^2) for a in (-1,1), and stable for every delay for
% % a <= -1. For r = 2.9 the minimum is attained at (m,n) = (1,1), a = 0.9:
% %   stable iff tau < tau* = 1.034726,
% % the threshold of the 1D version (r = 1.9, library example 43), which
% % Fridman, Orlov 2009, Remark 4, show to be tight [read].
% %
% % Certified (PIETOOLS_stability_sop, 'P' form, lpisettings_sop('light'),
% % MOSEK, 11/17/2026), tau = 0.5*tau*: infeasible with the stock light
% % settings (1,668,608 decision variables) and with the proposed light
% % row of the settings proposal of 11/15/2026 (1,680,004); no
% % certificate at this tier.
% %                                                                         % MMP, 11/17/2026
% % No second output st: no settings measured certify this example (gap     % MMP, 11/17/2026
% % study 11/17/2026, gap_results.md G6: infeasible at 2,514,256            % MMP, 11/17/2026
% % variables; the direct form exceeds 3,000,000 variables).                % MMP, 11/17/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: no settings measured certify the example, so
%                   no second output st (application-specific settings);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2 s3

%%% Executive Function:
evalin('base','stability = 1;');

r = 2.9;    tau = 0.9*1.034726;
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
x1 = pde_var([s1;s2],[0,pi;0,pi]);
x2 = pde_var([s1;s2;s3],[0,pi;0,pi;0,tau],[2;2;1]);   % 2nd order in s1, s2; 1st in s3
PDE_t = [diff(x1,'t')==diff(x1,s1,2)+diff(x1,s2,2)+r*x1-subs(x2,s3,tau);
         diff(x2,'t')==-diff(x2,s3);
         subs([x1;x2],s1,0)==0;     subs([x1;x2],s1,pi)==0;
         subs([x1;x2],s2,0)==0;     subs([x1;x2],s2,pi)==0;
         subs(x2,s3,0)==x1];


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
