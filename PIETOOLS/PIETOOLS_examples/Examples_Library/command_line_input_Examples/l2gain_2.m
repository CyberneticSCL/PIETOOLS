% L2GAIN_2 declares a second-order PDE on [0,1] with distributed
% disturbance s*w(t), BCs x(0) = 0, x_s(1) = 0, and output
% z(t) = int_0^1 x(t,s) ds; displays it, converts it to a PIE and displays
% the PIE. Workspace: pde ('pde_struct') and pie ('pie_struct').

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - l2gain_2
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
% MMP, 11/19/2026: Release comment pass (comment-only).

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
%   PDE: x_{t} = x_{ss} + s*w(t)                    |                       (gamma = 0.3333)
%   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |
%   Out: z(t) = int(x(t,s),s,0,1)                   |         
pde_var x output z input w
x.vars = s;     x.dom = [0,1];  
%% Define equations
eq_PDE = [diff(x,t)==-diff(x,s,2)+s*w]; %   PDE: x_{t} = x_{ss} + s*w(t)
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %   BCs: x(s=0) = 0, x_{s}(s=1) = 0
eq_out = [z==int(x,s,[0,1])]; %   Out: z(t) = int(x(t,s),s,0,1)
%% initialize pde system;
pde = initialize([eq_PDE;eq_BC;eq_out]);
%% display pde to verify and convert to pie
display(pde);
pie = convert(pde,'pie');
display(pie);



% %% Define dependent variables and system variable
% %   PDE: x_{t} = x_{ss} + s*w(t)                    |                       (gamma = 0.3333)
% %   BCs: x(s=0) = 0,        x_{s}(s=1) = 0          |
% %   Out: z(t) = int(x(t,s),s,0,1)                   |         
% x=state('pde'); z = state('out'); w=state('in');
% pde = sys();     
% %% Define equations
% eq_PDE = [diff(x,t)==-diff(x,s,2)+s*w]; %   PDE: x_{t} = x_{ss} + s*w(t)
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %   BCs: x(s=0) = 0, x_{s}(s=1) = 0
% eq_out = [z==int(x,s,[0,1])]; %   Out: z(t) = int(x(t,s),s,0,1)
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC;eq_out]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
