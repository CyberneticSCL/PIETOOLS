function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom(GUI,params)  % MMP, 11/19/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings (this example has none).
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         example-specific settings (lpisettings_sop structure)       % MMP, 11/19/2026
%               for the tests this file names                               % MMP, 11/19/2026
%
% %---------------------------------------------------------------------% %
% % 2D heat equation with a 1D distributed disturbance and a 2D output,
% % with exact L2 gain (Lidstrom, Rantzer, Morris 2016, Theorem 1 with
% % B = 0):
% % PDE         x_{t}(t,s1,s2) = x_{s1s1} + x_{s2s2} + w(t,s1);
% %                                               (s1,s2) in [0,1]^2
% % OUT         z(t,s1,s2) = x(t,s1,s2);
% % With BCs    x = 0 at s1 = 0, 1 and at s2 = 0, 1.
% %
% % Exact gain: A (the Dirichlet Laplacian) is self-adjoint and strictly
% % negative and C = I, so the gain is attained at frequency 0 and equals
% % ||A^{-1} H|| [read, Theorem 1]. Separation in sin(n*pi*s1) gives
% % gamma* = sup_n ||psi_n||_{L2(0,1)}, psi_n'' - (n*pi)^2*psi_n = 1,
% % psi_n(0) = psi_n(1) = 0, attained at n = 1 with
% % psi_1(s) = -(1/pi^2)*(1 - cosh(pi*(s-1/2))/cosh(pi/2)) [derived]:
% %   gamma* = 0.04571818.
% % Finite differences (5-point Laplacian, N interior nodes per direction):
% % 0.04571809 (N = 40), 0.04571814 (N = 60).
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier bounds the gain within 0.5 %,
%                   so no second output st (application-specific settings);
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
% function [PDE_t] = PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom(GUI,params)  % MMP, 11/19/2026 (was)

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_gain = 1;');
evalin('base','Hinf_gain_form = ''P'';');  % form of the gain test          % MMP, 11/19/2026

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
x = pde_var([s1;s2],[0,1;0,1]);
w = pde_var('in',1,s1,[0,1]);                   % 1D disturbance, constant in s2
z = pde_var('out',1,[s1;s2],[0,1;0,1]);
PDE_t = [diff(x,'t')==diff(x,s1,2)+diff(x,s2,2)+w;
         z==x;
         subs(x,s1,0)==0;           subs(x,s1,1)==0;
         subs(x,s2,0)==0;           subs(x,s2,1)==0];


% Example-specific settings                                                 % MMP, 11/19/2026
st = lpisettings_sop('light');                                              % MMP, 11/19/2026
st.lpi.pos.len_com{2} = [2 1];                                              % MMP, 11/19/2026
st.lpi.pos.len_sos{2} = [1 1;1 1;0 0];                                      % MMP, 11/19/2026
st.lpi.pos.patterns{2} = 'uniform';                                         % MMP, 11/19/2026
st.lpi.ineq.wcross = 1;                                                     % MMP, 11/19/2026

if GUI~=0
    disp('No GUI representation available for this system.')
end

end
% @inproceedings{lidstrom2016hinf,
%   title={H-infinity optimal control for infinite-dimensional systems with strictly negative generator},
%   author={Lidstr{\"o}m, Carolina and Rantzer, Anders and Morris, Kirsten},
%   booktitle={Proceedings of the 55th IEEE Conference on Decision and Control},
%   pages={5275--5280},
%   year={2016},
%   doi={10.1109/CDC.2016.7799077}
% }
