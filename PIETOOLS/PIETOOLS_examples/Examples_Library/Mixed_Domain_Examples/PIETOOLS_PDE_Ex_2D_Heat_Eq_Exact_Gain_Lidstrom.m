function [PDE_t,st] = PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom(GUI,params)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PIETOOLS PDE Examples
% [PDE_t,ST] = PIETOOLS_PDE_EX_2D_HEAT_EQ_EXACT_GAIN_LIDSTROM(GUI,PARAMS)
% returns the 2-D heat equation below, with a disturbance w(t,s1) on one
% variable and the full-state output z = x, as a 'pde_struct' object, and
% the recommended settings ST of its H-infinity gain test. It sets
% Hinf_gain = 1 and Hinf_gain_form = 'P' (storage <T x, R T x>,
% R >= eppos I) in the base workspace; a nonzero GUI prints that no GUI
% representation exists.
%
% INPUT
% - GUI:        Binary index {0,1} indicating whether or not a GUI
%               implementation of the example should be produced.
% - params:     Optional parameters for the example, that should be
%               specified as a cell of strings (this example has none).
%
% OUTPUT
% - PDE_t:      PDE structure defining the example system in the term-based
%               format.
% - st:         recommended settings (lpisettings structure) for the
%               test this file names.
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
% %---------------------------------------------------------------------% %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% Initial coding MMP - 11/17/2026

% Determine the location of this example file <-- DO NOT MOVE THE FILE
loc = mfilename('fullpath');
root = fileparts(loc);

% Initialize variables
pvar s1 s2

%%% Executive Function:
evalin('base','Hinf_gain = 1;');
evalin('base','Hinf_gain_form = ''P'';');  % form of the gain test

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


% Example-specific settings: the light tier, with the storage on the 2-D
% space from the basis operators aa and bb only ('uniform'), s-lengths
% [2 1] of the integral b and theta-lengths [1 1] of a and b; slack weight
% degree 1 in the directions a space does not own (wcross)
st = lpisettings('light');
st.lpi.pos.len_com{2} = [2 1];
st.lpi.pos.len_sos{2} = [1 1;1 1;0 0];
st.lpi.pos.patterns{2} = 'uniform';
st.lpi.ineq.wcross = 1;

% No GUI representation of this system
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
