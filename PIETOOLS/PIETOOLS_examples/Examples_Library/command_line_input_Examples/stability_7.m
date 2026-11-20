%   PDE:  x_{t} = a(s)*x_{ss}                    | a = s^3 - s^2 + 2     (unstable for lam > 4.66)                                                                                                      
%                 + b(s)*x_{s}                          | b = 3*s^2 - 2*s                 exponentially stable and FE PDE stable for lam < 4.66 
%                 + c(s,lam)*x                          | c =-0.5*s^3 + 1.3*s^2            Gahlawat 2017 [4]
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0     |    - 1.5*s + 0.7 +lam
%                                                              | lam = 4.65   
%
% STABILITY_7 declares the parabolic PDE above, with spatially varying
% coefficients a, b, c and lam = 4.65, x(0) = 0, x_s(1) = 0; displays it,
% converts it to a PIE and displays the PIE. Workspace: PDE ('pde_struct')
% and PIE ('pie_struct'). If a cell array params of assignment strings
% exists in the workspace, its entries are evaluated after the default
% parameters (c is computed from lam before that).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_7
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
% Define dependent variables and system variables
x = pde_var(s,[0,1]);
%% Specify the parameters
lam = 4.65; a = s^3-s^2+2; b=3*s^2-2*s; c=-0.5*s^3+1.3*s^2-1.5*s+0.7+lam;
 if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
 end
%% Define equations
eq_PDE = diff(x,t)==c*x+b*diff(x,s)+a*diff(x,s,2); % 	PDE: x_{t} = a(s)*x_{ss}+b(s)*x_{s}+c(s,lam)*x
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE:  x_{t} = a(s)*x_{ss}                       | a = s^3 - s^2 + 2     (unstable for lam > 4.66)            
% %                 + b(s)*x_{s}                      | b = 3*s^2 - 2*s               Gahlawat 2017 [4]
% %                 + c(s,lam)*x                      | c =-0.5*s^3 + 1.3*s^2 
% %   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            |    - 1.5*s + 0.7 +lam
% %                                                   | lam = 4.66
% x = state('pde');
% pde = sys(); lam = 4.66; a = s^3-s^2+2;b=3*s^2-2*s;c=-0.5*s^3+1.3*s^2-1.5*s+0.7+lam;
% %% Define equations
% eq_PDE = diff(x,t)==c*x+b*diff(x,s)+a*diff(x,s,2); % 	PDE: x_{t} = a(s)*x_{ss}+b(s)*x_{s}+c(s,lam)*x
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,  x_{s}(s=1) = 0
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
