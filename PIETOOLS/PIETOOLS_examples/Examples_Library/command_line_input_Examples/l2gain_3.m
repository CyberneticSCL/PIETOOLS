% L2GAIN_3 declares the parabolic PDE with spatially varying coefficients
%   x_t = a(s)*x_ss + b(s)*x_s + c(s,lam)*x + w(t),  s in [0,1],
%   x(0) = 0,   x_s(1) = 0,   z(t) = x(t,1),
% with a, b, c and lam = 4.66 as below; displays it, converts it to a PIE
% and displays the PIE. Workspace: pde ('pde_struct') and pie
% ('pie_struct').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - l2gain_3
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
% Initial coding SS - 11/17/2022
% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
%   PDE:  x_{t} = a(s)*x_{ss}                       | a = s^3 - s^2 + 2     (unstable for lam > 4.66)            
%                 + b(s)*x_{s}                      | b = 3*s^2 - 2*s               Gahlawat 2017 [4]
%                 + c(s,lam)*x +w(t)                | c =-0.5*s^3 + 1.3*s^2 
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            |    - 1.5*s + 0.7 +lam
%   Out: z(t) = x(t,1)                              | lam = 4.66
pde_var x output z input w
x.vars = s;     x.dom = [0,1]; 
lam = 4.66; a = s^3-s^2+2; b=3*s^2-2*s; c=-0.5*s^3+1.3*s^2-1.5*s+0.7+lam;
%% Define equations
eq_PDE = diff(x,t)==c*x+b*diff(x,s)+a*diff(x,s,2)+w; % 	PDE: x_{t} = a(s)*x_{ss}+b(s)*x_{s}+c(s,lam)*x+w
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
eq_out = [z==subs(x,s,1)]; % Out: z(t) = x(t,1)
%% initialize pde system;
pde = initialize([eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);



% %% Define dependent variables and system variable
% %   PDE:  x_{t} = a(s)*x_{ss}                       | a = s^3 - s^2 + 2     (unstable for lam > 4.66)            
% %                 + b(s)*x_{s}                      | b = 3*s^2 - 2*s               Gahlawat 2017 [4]
% %                 + c(s,lam)*x +w(t)                | c =-0.5*s^3 + 1.3*s^2 
% %   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            |    - 1.5*s + 0.7 +lam
% %   Out: z(t) = x(t,1)                              | lam = 4.66
% x = state('pde'); w = state('in'); z = state('out');
% pde = sys(); lam = 4.66; a = s^3-s^2+2;b=3*s^2-2*s;c=-0.5*s^3+1.3*s^2-1.5*s+0.7+lam;
% %% Define equations
% eq_PDE = diff(x,t)==c*x+b*diff(x,s)+a*diff(x,s,2)+w; % 	PDE: x_{t} = a(s)*x_{ss}+b(s)*x_{s}+c(s,lam)*x+w
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
% eq_out = [z==subs(x,s,1)]; % Out: z(t) = x(t,1)
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC; eq_out]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
