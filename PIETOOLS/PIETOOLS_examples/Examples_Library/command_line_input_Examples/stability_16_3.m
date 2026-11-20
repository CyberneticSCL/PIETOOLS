% STABILITY_16_3 declares the beam model, fourth order in time,
%   0 = alp*w_ssss + bet*w_tt - gam*w_ttss + delt*w_tttt,
% in the states x1 = w_t, x2 = w_ttt, x3 = w_tt, x4 = w on [0,1], with the
% six BCs of eq_BC below (one of them non-local) and alp = bet = delt = 1,
% gam = 2; displays it, converts it to a PIE and displays the PIE.
% Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell array
% params of assignment strings (e.g. {'gam=1;'}) exists in the workspace,
% its entries are evaluated after the default parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_16_3
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
% Initial coding DB - 06/19/2026

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%   PDE: 0 = alp*w_{ssss} + bet*w_{tt}              |
%             - gam*w_{ttss} + delt w_{tttt}             |
%     BCs  w(s=0) = 0,        w_{sss}(s=0)-w_{stt}(s=0)+w_s(s=0)=0, 
%                 w_{ss}(s=1)=beta*wtt(s=1), alp*w_{sss}(s=1)=gam*w_{stt}(s=1) 
% alp=E*I, beta= rho*A, gam= rhoI+E*I*rho/k/G, rho^2*I/k/G 
%% Specify the parameters
% Assuming E=I=rho=A=G=k=1,
alp = 1;     bet = 1;     gam = 2;     delt=1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%   Use states:                                     |
%        x1 = w_{t},    x2 = w_{ttt},              |
%        x3 = w_{tt},       x4 = w.                 |
%       =>                                       
%   PDE: x1_{t} = x3,
%             x2_{t} = -alp/delt*x4_{ssss} - bet/delt*x3           |
%                 + gam/delt*x3_{ss}                     |
%       ,  x3_{t} = x2,  x4_{t} = x1    |
%   BCs: x4(s=0) = 0,        x_4_{sss}-x3_{s}+x4_{s}(s=0) = 0          |
%        x4_{ss}(s=1) - x3(s=1) = 0                   |
%        alp*x4_{sss}(s=1) - gam*x3_{s}(s=1) = 0              |
%        x3(s=0) = 0    x3_{s}(s=0) = 0        |
%% Define dependent variables and system variable
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;   
%% Define equations
eq_PDE = [diff(x1,t,1)==x3; %  PDE:  x1_{t} = x3,
                    diff(x2,t,1)==-alp/delt*diff(x4,s,4)-bet/delt*x3+gam/delt*diff(x3,s,2); % x2_{t} = -alp/delt*x4_{ssss} - bet/delt*x3 + gam/delt*x3_{ss}  
                  diff(x3,t,1)==x2;%x3_{t} = x2,   
                  diff(x4,t,1)==x1;];% x4_{t} = x1    |
eq_BC = [subs(x4,s,0)==0; % BCs: x4(s=0) = 0,  
                 subs(diff(x4,s,3),s,0)- subs(diff(x3,s,1),s,0)+ subs(diff(x4,s,1),s,0)==0;%x_4_{sss}-x3_{s}+x4_{s}(s=0) = 0  
                subs(diff(x4,s,2),s,1)-subs(x3,s,1)==0; %   x4_{ss}(s=1) - x3(s=1) = 0                         |
                subs(diff(x4,s,3),s,1)- gam*subs(diff(x3,s,1),s,1)==0;   %  x4_{sss}(s=1) - gam*x3_{s}(s=1) = 0              |
                subs(x3,s,0)==0; %        x3(s=0) = 0         |
                 subs(diff(x4,s,1),s,0) + bet*int(x3,s,0,1)==0]; % non-local condition
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);

