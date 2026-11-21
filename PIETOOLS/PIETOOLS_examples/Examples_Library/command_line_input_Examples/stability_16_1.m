% STABILITY_16_1 declares a beam model, fourth order in time, in the states
% x1 = w_ttt, x2 = w_t, x3 = w_tt, x4 = w on [0,1]:
%   x1_t = -alp*x4_ssss - bet*x3 + gam*x3_ss,   x2_t = x3,   x3_t = x1,
%   x4_t = x2,
% with the BCs of eq_BC below (one of them non-local), alp = bet = 1,
% gam = 2; displays it, converts it to a PIE and displays the PIE.
% Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell array
% params of assignment strings (e.g. {'gam=1;'}) exists in the workspace,
% its entries are evaluated after the default parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_16_1
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
%
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi) 
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
% 
%% Specify the parameters
alp = 1;     bet = 1;     gam = 2;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   We rewrite the system as a single equation:     |                       4th order implementation (stable)
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}              |
%             - gam*w_{ttss} + w_{tttt}             |
%   Use states:                                     |
%        x1 = w_{ttt},      x2 = w_{t},             |
%        x3 = w_{tt},       x4 = w.                 |
%       =>                                        
%   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3           |
%                 + gam*x3_{ss}                     |
%        x2_{t} = x3,  x3_{t} = x1,  x4_{t} = x2    |
%   BCs: x4(s=0) = 0,        x4_{s}(s=0) = 0          |
%        x4_{ss}(s=1) - x4(s=1) = 0                   |
%        x4_{sss}(s=1) - x4_{s}(s=1) = 0              |
%        x2(s=0) = 0     x2_{s}(s=0) = 0         |
%        x3(s=0) = 0    x3_{s}(s=0) = 0        |
%% Define dependent variables and system variable
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;   
%% Define equations
eq_PDE = [diff(x1,t,1)==-alp*diff(x4,s,4)-bet*x3+gam*diff(x3,s,2); %   PDE: x1_{t} = -alp*x4_{ssss} - bet*x3   + gam*x3_{ss}         |
                  diff(x2,t,1)==x3; %    x2_{t} = x3,  
                  diff(x3,t,1)==x1;% x3_{t} = x1,  
                  diff(x4,t,1)==x2;];% x4_{t} = x2    |
% BCs: x1(0) = x2(0) = x3(0) = x4(0) = 0 and the four below
eq_BC = [ subs(x1,s,0)==0; subs(x2,s,0)==0;
                subs(x3,s,0)==0;subs(x4,s,0)==0; 
                subs(diff(x4,s,1),s,0)+int(x3,s,0,1)==0;% x4_{s}(s=0) + int_0^1 x3 ds = 0
                subs(diff(x4,s,3),s,0)-subs(diff(x3,s,1),s,0)+subs(diff(x4,s,1),s,0)==0; % x4_{sss}(s=0) - x3_{s}(s=0) + x4_{s}(s=0) = 0
                subs(diff(x4,s,2),s,1)- subs(x3,s,1)==0;   % x4_{ss}(s=1) - x3(s=1) = 0
                subs(diff(x4,s,3),s,1)-gam*subs(diff(x3,s,1),s,1)==0]; % x4_{sss}(s=1) - gam*x3_{s}(s=1) = 0
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
