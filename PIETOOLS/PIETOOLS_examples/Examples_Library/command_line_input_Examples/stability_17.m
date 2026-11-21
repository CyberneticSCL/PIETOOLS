% STABILITY_17 declares the tip-damped wave equation in the states
% x1 = u_s, x2 = u_t on [0,1]:
%   x1_t = x2_s,   x2_t = x1_s,   x2(0) = 0,   x1(1) + k*x2(1) = 0,
% with k = 1; displays it, converts it to a PIE and displays the PIE.
% Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell array
% params of assignment strings (e.g. {'k=2;'}) exists in the workspace,
% its entries are evaluated after the default parameters.
%
% Tip-damped wave equation
%   PDE: u_{tt} = u_{ss}                            | k=1    
%   BCs: u(s=0) = 0,                                |               Peet 2019 [8] (Example 8.2)
%        u_{s}(s=1) = -k*u_{t}(s=1)                 |
%   Use states x1 = u_{s}, x2 = u_{t}.              |     (PDE exponentially stable for k>0)                 
%       =>                                          |                      (Finite Energy PDE stable  for k>0) 
%   PDE: x1_{t} = x2_{s},   x2_{t} = x1_{s}         |
%   BCs: x2(0) = 0,         x1(1) + k*x2(1) = 0     |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_17
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
clc; clearvars -except params; clear stateNameGenerator
pvar s t; % define independent variables
%% Specify the parameters
k = 1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable
pde_var x1 x2
x1.vars = s;    x2.vars = s; 
%% Define equations
eq_PDE = [diff(x1,t,1)==diff(x2,s,1); 
                  diff(x2,t,1)==diff(x1,s,1);];
eq_BC = [subs(x2,s,0)==0;subs(x1,s,1)+k*subs(x2,s,1)==0;];
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);