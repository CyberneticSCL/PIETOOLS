% L2GAIN_23 declares the parabolic PDE with spatially varying coefficients
%   x_t = A2(s)*x_ss + A1(s)*x_s + A0(s)*x + w(t),  s in [0,1],
%   x(0) = 0,   x_s(1) = 0,   z(t) = x(t,1),
% with A2, A1, A0 and lam = 4.6 as below; displays it, converts it to a
% PIE and displays the PIE. Workspace: PDE ('pde_struct') and PIE
% ('pie_struct').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - l2gain_23
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
% Initial coding DB - 05/31/2026

%clc; clear; clear stateNameGenerator
pvar s t; % define independent variables

%% Define dependent variables and system variable
%   PDE: x_{t} = A2(s)*x_{ss} + A1(s)*x_{s} + A0(s)*x + w(t)
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0
%   Out: z(t) = x(t,1)
%   A2 = s^3 - s^2 + 2;   A1 = 3*s^2 - 2*s
%   A0 = -0.5*s^3 + 1.3*s^2 - 1.5*s + 0.7 + lam;  lam = 4.6
pde_var x output z input w
x.vars = s;     x.dom = [0,1];
lam = 4.6;
A2 = s^3-s^2+2;
A1 = 3*s^2-2*s;
A0 = -0.5*s^3+1.3*s^2-1.5*s+0.7+lam;

%% Define equations
eq_PDE = diff(x,t)==A0*x+A1*diff(x,s)+A2*diff(x,s,2)+w;
eq_BC = [subs(x,s,0)==0; subs(diff(x,s),s,1)==0];
eq_out = z==subs(x,s,1);

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC; eq_out]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
