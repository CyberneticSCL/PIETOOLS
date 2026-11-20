% STABILITY_15 declares the Timoshenko beam in the hyperbolic/diffusive
% states x1 = w_t, x2 = w_s, x3 = phi_t, x4 = phi on [0,1]:
%   x1_t = x2_s - x4_s,   x2_t = x1_s,   x3_t = x2 - x4 + x4_ss,   x4_t = x3,
%   x4(0) = x3(0) = x1(0) = 0,   x4_s(1) = 0,   x2(1) - x4(1) = 0,
% with all coefficients 1 and no damping (c = 0); displays it, converts it
% to a PIE and displays the PIE. Workspace: PDE ('pde_struct') and PIE
% ('pie_struct'). The defaults k, aa, II, g are not used by the equations.
% If a cell array params of assignment strings exists in the workspace,
% its entries are evaluated after the defaults.
%
% Timoshenko Beam
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi) 
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
%   Assume all parameters are 1, and use states:    | Hyperbolic/diffusive implementation (c = 0: stable, not exponentially stable)
%        x1 = w_{t},    x2 = w_{s},                 |
%        x3 = phi_{t},  x4 = phi.                   |                  (PDE stable with damping c >0)
%       =>                                          |
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
%        x3_{t} = x2 - x4 -c x3                         |                   % MMP, 11/19/2026 (was)
%        x3_{t} = x2 - x4 + x4_{ss} - c x3              |                   % MMP, 11/19/2026
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_15
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
% MMP, 11/19/2026: the bending term x4_{ss} (E*II*phi_{ss}) restored in
%                   x3_t, as in PIETOOLS_PDE_Ex_Timoshenko_Beam_2.m (Peet
%                   et al. 2019, Sec. 8.1); at c = 0 the beam is stable, not
%                   exponentially stable (Peet, Automatica 2021, Sec. 11.1).
% MMP, 11/19/2026: Release comment pass (comment-only).

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Specify the parameters 
k = 1;    aa = 1;    II = 1;    g = 1;   
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable
% use PDE states:  (Hyperbolic/diffusive implementation)
%        x1 = w_{t},    x2 = w_{s},                 |
%        x3 = phi_{t},  x4 = phi.                   |
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
% %        x3_{t} = x2 - x4                           |                     % MMP, 11/19/2026 (was)
%        x3_{t} = x2 - x4 + x4_{ss}                 |                       % MMP, 11/19/2026
%        x4_{t} = x3                                |
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0           
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;   
%% Define equations
%   PDE: x1_{t} = x2_{s} - x4_{s}                   |
%        x2_{t} = x1_{s}                            |
% %        x3_{t} = x2 - x4                           |                     % MMP, 11/19/2026 (was)
%        x3_{t} = x2 - x4 + x4_{ss}                 |                       % MMP, 11/19/2026
%        x4_{t} = x3                                |
eq_PDE = [diff(x1,t)==diff(x2,s)-diff(x4,s); 
          diff(x2,t)==diff(x1,s);
%         diff(x3,t)==x2-x4;                                                % MMP, 11/19/2026 (was)
          diff(x3,t)==x2-x4+diff(x4,s,2);                                   % MMP, 11/19/2026
          diff(x4,t)==x3];
%   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
%        x3(0) = 0,     x1(0) = 0,                  |
%        x2(1) - x4(1) = 0 
eq_BC = [subs(x4,s,0)==0;subs(x3,s,0)==0;subs(diff(x4,s),s,1)==0;
         subs(x1,s,0)==0;subs(x2,s,1)-subs(x4,s,1)==0]; 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: x1_{t} = x2_{s} - x4_{s}                   |
% %        x2_{t} = x1_{s}                            |
% %        x3_{t} = x2 - x4                           |
% %        x4_{t} = x3                                |
% %   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
% %        x3(0) = 0,     x1(0) = 0,                  |
% %        x2(1) - x4(1) = 0           
% x1=state('pde');x2=state('pde');x3=state('pde');x4=state('pde'); 
% pde = sys();     
% %% Define equations
% %   PDE: x1_{t} = x2_{s} - x4_{s}                   |
% %        x2_{t} = x1_{s}                            |
% %        x3_{t} = x2 - x4                           |
% %        x4_{t} = x3                                |
% eq_PDE = [diff(x1,t)==diff(x2,s)-diff(x4,s); 
%           diff(x2,t)==diff(x1,s);
%           diff(x3,t)==x2-x4;
%           diff(x4,t)==x3];
% %   BCs: x4(0) = 0,     x4_{s}(1) = 0,              |
% %        x3(0) = 0,     x1(0) = 0,                  |
% %        x2(1) - x4(1) = 0 
% eq_BC = [subs(x4,s,0)==0;subs(x3,s,0)==0;subs(diff(x4,s),s,1)==0;
%          subs(x1,s,0)==0;subs(x2,s,1)-subs(x4,s,1)==0]; 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
