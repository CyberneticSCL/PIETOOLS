% STABILITY_10 declares the heat equation coupled with an ODE at the
% boundary,
%   xo_t = xo + x_s(0),   x_t = x_ss on [0,1],   x(0) = -xo,   x(1) = k*xo,
% with k = -2; displays it, converts it to a PIE and displays the PIE.
% Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell array
% params of assignment strings (e.g. {'k=-3;'}) exists in the workspace,
% its entries are evaluated after the default parameters.
%
%   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2                (stable for k=-2)            
%   PDE: x_{t} = x_{ss}                             |                               Tang 2011 [13]
%   BCs: x(s=0) = xo,    x(s=1) = k*xo              |                           (PDE stable )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_10
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
% Define dependent variables and system variable
x=pde_var(); X = pde_var(s,[0,1]);
%% Specify the parameters
k = -2;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_ODE=diff(x,t)==x+subs(diff(X,s),s,0); % ODE: xo_{t} = xo + x_{s}(s=0)
eq_PDE = diff(X,t)==diff(X,s,2); % 	PDE: x_{t} = x_{ss}
eq_BC = [subs(X,s,0)==-x;subs(X,s,1)==k*x]; %  BCs: x(s=0) = -xo,    x(s=1) = k*xo
%% initialize pde system;
PDE = initialize([eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   ODE: xo_{t} = xo + x_{s}(s=0)               	| k = -2   
% %   PDE: x_{t} = x_{ss}                             |                            
% %   BCs: x(s=0) = -xo,    x(s=1) = k*xo             |
% x=state('ode'); X = state('pde');
% pde = sys(); 
% k = -2;
% %% Define equations
% eq_ODE=diff(x,t)==x+subs(diff(X,s),s,0); % ODE: xo_{t} = xo + x_{s}(s=0)
% eq_PDE = diff(X,t)==diff(X,s,2); % 	PDE: x_{t} = x_{ss}
% eq_BC = [subs(X,s,0)==-x;subs(X,s,1)==k*x]; %  BCs: x(s=0) = -xo,    x(s=1) = k*xo
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_ODE;eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
