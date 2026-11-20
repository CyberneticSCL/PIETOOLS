% STABILITY_13_3 declares the Euler-Bernoulli beam with damping in the
% states x1 = w_t, x2 = w on [0,1]:
%   x1_t = -alp/bet*x2_ssss - c*x1,   x2_t = x1,
%   x2(0) = x2_s(0) = 0,   x2_ss(1) = x2_sss(1) = 0,
% with alp = 0.1, bet = 1, c = 0; displays it, converts it to a PIE and
% displays the PIE. Workspace: PDE ('pde_struct') and PIE2 ('pie_struct').
% If a cell array params of assignment strings (e.g. {'c=0.1;'}) exists in
% the workspace, its entries are evaluated after the default parameters.
%
%% Euler-Bernoulli Equation
% using PDE states x1=w_{t}, x2=w and PIE states (w_{t},w_{ssss})
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%   We rewrite the system as a single equation:     |                     
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt} + c*bet*w_{t} |
%     BCs w(s=0) = 0,        w_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
%        w_{s}(s=0) = 0,    w_{sss}(s=1) = 0        |
%   Use states:                                     |
%        x1 = w_{t},      x2 = w.                 |
%       =>                                       
%   PDE: x1_{t} = -alp/bet*x2_{ssss},
%             x2_{t} = x_1           |
%   BCs: x2(s=0) = 0,        x2_{ss}(s=1) = 0      |
%        x2_{s}(s=0)= 0        x2_{sss}(s=1)= 0    |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_13_3
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
%clear stateNameGenerator
pvar s t; % define independent variables
% alp = E*I, bet = rho*A
%% Specify the parameters
% Assuming E=I=rho=A=G=k=1,
alp = .1;     bet = 1;
% PDE states x1, x2; damping coefficient c
pde_var x1 x2
c=0;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
x1.vars = s;   x2.vars = s;
%% Define equations
eq_PDE = [diff(x1,t,1)==-alp/bet*diff(x2,s,4)-c*x1; %  %   PDE: x1_{t} = -alp/bet*x2_{ssss}-c*x1,
                    diff(x2,t,1)==x1;]; %             x2_{t} = x_1           |
eq_BC = [subs(x2,s,0)==0;
    subs(diff(x2,s,2),s,1)==0;%%   BCs: x2(s=0) = 0,        x2_{ss}(s=1) = 0      |
                subs(diff(x2,s,1),s,0)==0; 
                subs(diff(x2,s,3),s,1)==0;];%        x2_{s}(s=0)= 0        x2_{sss}(s=1)= 0    |
% Initialize the PDE
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE2 = convert(PDE,'pie');
display(PIE2);
