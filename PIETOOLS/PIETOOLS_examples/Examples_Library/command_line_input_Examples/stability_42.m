% STABILITY_42 declares an ODE coupled with a heat equation with delays,
% the delays modeled by the transport states x3 on [1,1+tau] and x4 on
% [0,1]x[1,1+tau] (below), with the BCs
%   x2_{s1}(0) = 0,  x2(1) = 0,  x3(1+tau) = x1,
%   x4_{s1}(0,s2) = 0,  x4(1,s2) = 0,  x4(s1,1+tau) = x2(s1);
% displays it, converts it to a PIE and displays the PIE. Workspace: PDE
% ('pde_struct') and PIE ('pie_struct').
%
%   ODE: x1_{t} = A*x1 + A1*x3(t,1) + B*x2(t,s1=0)
%   PDE: x2_{t} = x2_{s1s1} + a*x2 + a2*x4(t,s1,s2=1)
%        x3_{t} = x3_{s2}
%        x4_{t} = x4_{s2}
%   tau = 1; A = -1; A1 = 0.5; B = -1; a = 1; a2 = 0.5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_42
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
pvar s1 s2 t; % define independent variables
%% Specify parameters
tau = 1;
A = -1;
A1 = 0.5;
B = -1;
a = 1;
a2 = 0.5;
%% Define dependent variables and system variable
pde_var x1 x2 x3 x4
x2.vars = s1;       x2.dom = [0,1];
x3.vars = s2;       x3.dom = [1,1+tau];
x4.vars = [s1;s2];  x4.dom = [0,1;1,1+tau];

%% Define equations
eq_PDE = [diff(x1,'t')==A*x1+B*subs(x2,s1,0)+A1*subs(x3,s2,1);
          diff(x2,'t')==diff(x2,s1,2)+a*x2+a2*subs(x4,s2,1);
          diff(x3,'t')==diff(x3,s2);
          diff(x4,'t')==diff(x4,s2)];
eq_BC = [subs(diff(x2,s1),s1,0)==0;
         subs(x2,s1,1)==0;
         subs(x3,s2,1+tau)==x1;
         subs(diff(x4,s1),s1,0)==0;
         subs(x4,s1,1)==0;
         subs(x4,s2,1+tau)==x2];

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
