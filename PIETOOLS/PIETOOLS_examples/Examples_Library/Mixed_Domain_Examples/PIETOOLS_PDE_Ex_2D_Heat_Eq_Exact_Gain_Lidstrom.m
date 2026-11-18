function [PDE_t] = PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom(GUI,params)
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
% %
% % Certified upper bounds (PIETOOLS_Hinf_gain_sop, 'P' form, MOSEK,
% % 11/17/2026), optimal exits:
% %   - stock light:  0.0459440156 (+0.494 %, 209,445 decision variables);
% %   - proposed light row of the settings proposal of 11/15/2026:
% %                   0.0459459178 (+0.498 %, 210,006);
% %   - stock heavy:  0.0459431727 (+0.492 %, 837,690).
% %                                                                         % MMP, 11/17/2026
% % No second output st: the light tier gives the bounds above, within      % MMP, 11/17/2026
% % 0.5 % of the exact gain; the light tier of 11/17/2026 (wcross 0, 2-D    % MMP, 11/17/2026
% % faces at offset 1) was not rerun.                                       % MMP, 11/17/2026
% %---------------------------------------------------------------------% %
%
% MMP, 11/17/2026: Initial coding;
% MMP, 11/17/2026: Header: the light tier bounds the gain within 0.5 %,
%                   so no second output st (application-specific settings);

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_gain = 1;');

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
