function [PDE_t, st] = PIETOOLS_PDE_Ex_Transport_Edge_Diffusion(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings e.g. {'r2=5;'}.
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         (optional) settings structure that certifies the example
%               (application-specific settings, settings report
%               11/16/2026 Sec. 5.8).
%
% %---------------------------------------------------------------------% %
% % Transport on an edge driving a diffusion in the normal direction
% % (mixed-domain test suite, system S17_transport_edge_diffusion):
% % PDE         x1_{t}(t,s1)    = -x1_{s1};               s1 in [0,1]
% %             x2_{t}(t,s1,s2) = x2_{s2s2} + r2*x2;      (s1,s2) in [0,1]^2
% % With BCs    x1(t,s1=0) = 0;
% %             x2(t,s1,s2=0) = x1(t,s1);   x2(t,s1,s2=1) = 0.
% % x2 carries no derivative in s1 (order 0 in s1, 2 in s2).
% %
% % Parameter r2 can be set; default r2 = 0.9*r2*.
% %
% % Threshold [derived]: x1(t,.) = 0 for t >= 1 (transport with zero
% % inflow); from then on each slice x2(t,s1,.) obeys the Dirichlet heat
% % equation on [0,1] with reaction r2, with eigenvalues r2 - k^2*pi^2,
% % k >= 1:
% %   exponentially stable iff r2 < r2* = pi^2.
% %
% % Storage form 'P' (P >= eppos*I): no coercive PI operator P certifies    % MMP, 11/19/2026
% % this system, for any r2 [argued, not checked by a second reader].       % MMP, 11/19/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding; the system and the measurements of the
%                   mixed-domain suite (S17), with the settings that
%                   certify it as the second output (maintainer:
%                   application-specific settings where the tiers do not);
% MMP, 11/18/2026: st: Q degree through st.sop.deg_Q (the direct form no
%                   longer reads lpi.var.deg) and the slack lift degree
%                   at its minimum, 10% fewer variables (header);
% MMP, 11/19/2026: st updated to the per-example search of 11/19/2026
%                   (maintainer: no default settings): st.lpi.ineq.wcross = 1;
%                   the form of its tests named in the workspace.
% MMP, 11/19/2026: the settings notes of the header (11/17-11/19) removed: they
%                   recorded tests and files the user does not have, not the problem;
% MMP, 11/19/2026: the remaining settings notes removed: they recorded tests,
%                   not the problem;
% MMP, 11/19/2026: the dated test records of the header removed: they
%                   recorded tests and files the user does not have, not
%                   the problem; the 'P' obstruction kept, reworded;

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','stability = 1;');
evalin('base','stability_form = ''Q'';');  % form of the stability test     % MMP, 11/19/2026

r2 = 0.9*pi^2;
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
x1 = pde_var(1,s1,[0,1]);
x2 = pde_var(1,[s1;s2],[0,1;0,1],[0;2]);       % order 0 in s1, 2 in s2
PDE_t = [diff(x1,'t')==-diff(x1,s1);
         diff(x2,'t')==diff(x2,s2,2)+r2*x2;
         subs(x1,s1,0)==0;
         subs(x2,s2,0)==x1;         subs(x2,s2,1)==0];


% Application-specific settings (header): the light tier with Q of the      % MMP, 11/18/2026 (was)
% direct form at degree 2 in two variables; built only when requested       % MMP, 11/18/2026 (was)
% Example-specific settings                                                 % MMP, 11/19/2026
if nargout>1
    st = lpisettings_sop('light');
    st.lpi.ineq.dD = 0;                                                     % MMP, 11/18/2026
    st.sop.deg_Q{2} = 2;                                                    % MMP, 11/18/2026
%    st.lpi.var.deg{2} = 2;                                                 % MMP, 11/18/2026 (was)
    st.lpi.ineq.wcross = 1;                                                 % MMP, 11/19/2026
end

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
