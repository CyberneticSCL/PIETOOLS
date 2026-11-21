% STABILITY_13_2 declares the Euler-Bernoulli beam as one equation, second
% order in time, in the state x = w on [0,1]:
%   x_tt = -alp/bet*x_ssss,   x(0) = x_s(0) = 0,   x_ss(1) = x_sss(1) = 0,
% with alp = 0.1, bet = 1; displays it, converts it to a PIE (which adds
% the state x_t) and displays the PIE. Workspace: PDE ('pde_struct') and
% PIE ('pie_struct'). If a cell array params of assignment strings (e.g.
% {'alp=1;'}) exists in the workspace, its entries are evaluated after the
% default parameters.
%
%% Euler-Bernoulli Equation
%   We rewrite the system as a single equation:     |                     
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}      |
%     BCs w(s=0) = 0,        w_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
%        w_{s}(s=0) = 0,    w_{sss}(s=1) = 0        |
%   Use states:                                     |
%        x1 = w_{t},      x2 = w.                 |
%       =>                                       
%   PDE: x1_{t} = -alp/bet*x2_{ssss},
%             x2_{t} = x_1           |
%   BCs: x2(s=0) = 0,        x2_{ss}(s=1) = 0      |
%        x2_{s}(s=0)= 0        x2_{sss}(s=1)= 0    |
% using PDE states x1=w_{t}, x2=w and PIE states (w_{t,ssss},w_{ssss})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_13_2
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
% alp = E*I, bet = rho*A
%% Specify the parameters
% Assuming E=I=rho=A=G=k=1,
alp = 0.1;     bet = 1; 
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
% PDE state x = w on [0,1]
pde_var x
x.vars = s;
%% Define equations
eq_PDE = diff(x,t,2)==-alp/bet*diff(x,s,4); % PDE: x_{tt} = -alp/bet*x_{ssss}
eq_BC = [subs(x,s,0)==0; subs(diff(x,s,2),s,1)==0;% BCs: x(s=0) = 0,        x_{ss}(s=1) = 0
                subs(diff(x,s,1),s,0)==0; subs(diff(x,s,3),s,1)==0;];%      x_{s}(s=0) = 0,      x_{sss}(s=1) = 0
% Initialize the PDE
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
