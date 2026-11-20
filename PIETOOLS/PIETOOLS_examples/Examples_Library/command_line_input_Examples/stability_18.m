% STABILITY_18 declares the tip-damped wave equation with reaction term and
% viscous damping in the states x1 = u_t, x2 = u on [0,1]:
%   x1_t = -2*ad*x1 - ad^2*x2 + x2_ss,   x2_t = x1,
%   x1(0) = x2(0) = 0,   k*x1(1) + x2_s(1) = 0,
% with k = ad = 1; displays it, converts it to a PIE and displays the PIE.
% Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell array
% params of assignment strings (e.g. {'k=2;'}) exists in the workspace,
% its entries are evaluated after the default parameters.
%
% Tip-damped wave eq with reaction term and viscous damping
%     | PDE: u_{tt} + 2*ad*u_{t} = -ad^2*u + u_{ss} | k = 1     Datko 1986 [9] (Test 7.5c)
%     | BCs: u(s=0) = 0                             |
%     |      u_{s}(s=1) = -k*u_{t}(s=1)             |
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   Use states x1 = u_{t}, x2 = u.                  | k = 1         (PDE exponentially dual stable )                 
%       =>                                          | ad = 1                    (Finite Energy PIE to PDE stable )   
%   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
%   BCs: x1(0) = 0,     x2(0) = 0,                  |
%        k*x1(1) + x2_{s}(1) = 0                    |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_18
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

%clc; clear;
clear stateNameGenerator
pvar s t; % define independent variables
%% Specify parameters
k = 1;    ad = 1;
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
eq_PDE = [diff(x1,t,1)==-2*ad*x1-ad^2*x2+diff(x2,s,2); %   PDE: x1_{t} = -2*ad*x1 - ad^2*x2 + x2_{ss}      |
                  diff(x2,t,1)==x1;];% x2_{t} = x1 
eq_BC = [subs(x1,s,0)==0;subs(x2,s,0)==0;%x1(0) = 0,     x2(0) = 0,                  |
                k*subs(x1,s,1)+subs(diff(x2,s,1),s,1)==0;]; %        k*x1(1) + x2_{s}(1) = 0      
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);