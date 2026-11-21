% L2GAIN_4 declares the 2-component reaction-diffusion PDE
%   x_t = Cm*x + (1/R)*x_ss + s*[1;1]*w(t),  s in [0,1],
%   x(0) = 0,   x_s(1) = 0,   z(t) = int_0^1 x(t,s) ds,
% with Cm and R = 2.7 as below; displays it, converts it to a PIE and
% displays the PIE. Workspace: pde ('pde_struct') and pie ('pie_struct').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - l2gain_4
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
% 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss} + s*w(t)      | R = 2.7
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [1, 1.5; 5, 0.2]     
%   Out: z(t) = int(x(t,s),s,0,1)                  | 
pde_var x output z input w
x.vars = s;     x.dom = [0,1];  x.size = 2;
z.size = 2;
R = 2.7;Cm = [1, 1.5; 5, 0.2];
%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2)+s*[1;1]*w; % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
eq_out = [z==int(x,s,[0,1])];    %   Out: z(t) = int(x(t,s),s,0,1) 
%% initialize pde system;
pde = initialize([eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);



% %% Define dependent variables and system variable
% % 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss} + s*w(t)      | R = 2.7
% %   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [1, 1.5; 5, 0.2]     
% %   Out: z(t) = int(x(t,s),s,0,1)                  | 
% x=state('pde',2); w = state('in'); z = state('out',2);
% pde = sys(); 
% R = 2.7;Cm = [1, 1.5; 5, 0.2];
% %% Define equations
% eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2)+s*[1;1]*w; % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
% eq_out = [z==int(x,s,[0,1])];    %   Out: z(t) = int(x(t,s),s,0,1) 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
