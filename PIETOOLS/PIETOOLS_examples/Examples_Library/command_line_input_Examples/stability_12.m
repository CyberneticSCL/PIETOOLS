% STABILITY_12 declares the 3-component reaction-diffusion PDE
%   x_t = Cm(s)*x + (1/R)*x_ss on [0,1],   x(0) = 0,   x_s(1) = 0,
% with Cm = [0 0 0; s 0 0; s -s^2 0] and R = 30; displays it, converts it
% to a PIE and displays the PIE. Workspace: PDE ('pde_struct') and PIE
% ('pie_struct'). If a cell array params of assignment strings (e.g.
% {'R=21;'}) exists in the workspace, its entries are evaluated after the
% default parameters.
%
%   PDE:  x_{t} = Cm*x + (1/R)*x_{ss}               | R = 30            (stable for R<21 (not tight))
%   BCs:  x(s=0) = 0,     x_{s}(s=1) = 0            | Cm = [0 0 0;                  Ahmadi 2015 [5] Adapted Example B
%                                                   |       s 0 0;
%                                                   |       -s^2 0 0]
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_12
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
% MMP, 11/20/2026: clearvars -except params replaces clear, which deleted
%   params before the params loop read it.

% Clear the workspace
%clc; clear; clear stateNameGenerator                                       % MMP, 11/20/2026 (was)
clc; clearvars -except params; clear stateNameGenerator                     % MMP, 11/20/2026
pvar s t; % define independent variables
% Define dependent variables and system variable
x=pde_var(3,s,[0,1]); 
%% Specify the parameters
%Cm = [0,0,0;s,0,0;-s^2,0,0];     R = 21;
Cm = [0 0 0; s 0 0; s -s^2 0];      R = 30;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE= convert(PDE,'pie');
display(PIE);



% x=state('pde',3); 
% pde = sys(); 
% R = 21;Cm = [0,0,0;s,0,0;-s^2,0,0];
% %% Define equations
% eq_PDE = diff(x,t)==Cm*x+(1/R)*diff(x,s,2); % 	PDE: x_{t} = Cm*x + (1/R)*x_{ss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,1)==0]; %  BCs: x(s=0) = 0,     x_{s}(s=1) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
