% STABILITY_11 declares the 2-component reaction-diffusion PDE
%   x_t = Cm*x + (1/R)*x_ss on [0,1],   x(0) = x(1) = 0,
% with Cm and R as assigned below; displays it, converts it to a PIE and
% displays the PIE. Workspace: PDE ('pde_struct')
% and PIE ('pie_struct'). If a cell array params of assignment strings
% (e.g. {'R=2.5;'}) exists in the workspace, its entries are evaluated
% after the default parameters.
%
% 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.9               (stable for R<2.7 (not tight))
%   BCs:  x(s=0) = 0,     x(s=1) = 0                | Cm = [1, 1.5; 5, 0.2]         Ahmadi 2014 [6] Example D
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_11
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
x=pde_var(2,s,[0,1]); 
%% Specify the parameters
R = 2,9;Cm = [1, 1.5; 5, 0.2];
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
eq_BC = [subs(x,s,0)==0;subs(x,s,1)==0]; %  BCs: x(s=0) = 0,     x(s=1) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% % 	PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 2.7
% %   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [1, 1.5; 5, 0.2]     
% x=state('pde',2); 
% pde = sys(); 
% R = 2.7;Cm = [1, 1.5; 5, 0.2];
% %% Define equations
% eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
