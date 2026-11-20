% L2GAIN_26 declares the 3-component reaction-diffusion PDE
%   x_t = Cm(s)*x + (1/R)*x_ss + [s;s;s]*w(t),  s in [0,1],
%   x(0) = x(1) = 0,   z(t) = int_0^1 x(t,s) ds,
% with Cm and R = 21 - 1e-3 as below; displays it, converts it to a PIE
% and displays the PIE. Workspace: PDE ('pde_struct') and PIE
% ('pie_struct').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - l2gain_26
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

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables

%% Define dependent variables and system variable
%   PDE: x_{t} = Cm(s)*x + (1/R)*x_{ss} + [s;s;s]*w(t)
%   BCs: x(s=0) = 0,        x(s=1) = 0
%   Out: z(t) = int(x(t,s),s,0,1)
%   Cm = [0,0,0; s,0,0; s^2,-s^3,0];   R = 21 - 1e-3
pde_var x output z input w
x.vars = s;     x.dom = [0,1];     x.size = 3;
z.size = 3;
Cm = [0,0,0; s,0,0; s^2,-s^3,0];
R = 21-1e-3;

%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2)+s*ones(3,1)*w;
eq_BC = [subs(x,s,0)==0; subs(x,s,1)==0];
eq_out = z==int(x,s,[0,1]);

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC; eq_out]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
