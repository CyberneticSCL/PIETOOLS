% 	ODE:  xo_{t} = A * xo + Bxr * x_{s}(s=a)        | a = 0    b = 1        (stable only for lam < -5.7799, so unstable at lam = pi^2-1)
%   PDE:  x_{t}  = lam * x + x_{ss} + Bpv * xo      | lam = pi^2-1                  Das 2018 [7] (Example 2)
%   BCs:  x(s=a) = 0,     x(s=b) = 0                | A, Bxr, and Bpv fixed
%
% STABILITY_8 declares the heat equation coupled with a 4-state ODE above
% (2-component PDE state X on [0,1]), with lam = pi^2-1 and the fixed
% matrices A, Bxr, Bpv; displays it, converts it to a PIE and displays the
% PIE. Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell
% array params of assignment strings (e.g. {'lam=5;'}) exists in the
% workspace, its entries are evaluated after the default parameters.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_8
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
%   params before it was read; line 1 gives the stability range as coded.

% Clear the workspace
%clc; clear; clear stateNameGenerator                                       % MMP, 11/20/2026 (was)
clc; clearvars -except params; clear stateNameGenerator                     % MMP, 11/20/2026
pvar s t; % define independent variables
% Define dependent variables and system variable
x=pde_var(4); X = pde_var(2,s,[0,1]);
%% Secify the parameters
lam = pi^2-1; a = 0;b=1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
% Fixed ODE matrix A and coupling matrices Bpv (ODE to PDE), Bxr (PDE to ODE)
A = [-1.2142,  1.9649,  0.2232,  0.5616;
            -1.8042, -0.7260, -0.3479,  5.4355;
            -0.2898,  0.7381, -1.7606,  0.8294;
            -0.9417, -5.3399, -1.0704, -0.7590];
Bpv = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
Bxr = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
%% Define equations
eq_ODE=diff(x,t)==A*x+Bxr*subs(diff(X,s),s,a); % ODE: xo_{t} = A * xo + Bxr * x_{s}(s=0)
eq_PDE = diff(X,t)==lam*X+diff(X,s,2)+Bpv*x;   % PDE: x_{t}  = lam * x + x_{ss} + Bpv * xo
eq_BC = [subs(X,s,a)==0;subs(X,s,b)==0];       % BCs: x(s=a) = 0,     x(s=b) = 0 
%% initialize pde system;
PDE = initialize([eq_ODE;eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% % 	ODE:  xo_{t} = A * xo + Bxr * x_{s}(s=a)        | a = 0    b = 1        
% %   PDE:  x_{t}  = lam * x + x_{ss} + Bpv * xo      | lam = pi^2-1          
% %   BCs:  x(s=a) = 0,     x(s=b) = 0                | A, Bxr, and Bpv fixed
% x=state('ode',4); X = state('pde',2);
% pde = sys(); lam = pi^2-1; a = 1;b=1;
% A = [-1.2142,  1.9649,  0.2232,  0.5616;
%             -1.8042, -0.7260, -0.3479,  5.4355;
%             -0.2898,  0.7381, -1.7606,  0.8294;
%             -0.9417, -5.3399, -1.0704, -0.7590];
% Bpv = [-2.5575 0 1.0368 0;-1.8067 0.4630 1.3621 0];
% Bxr = [-1.5368 0;0 0.8871;1.0656 0;1.1882 0];
% %% Define equations
% eq_ODE=diff(x,t)==A*x+Bxr*diff(subs(X,s,0),s); % ODE: xo_{t} = A * xo + Bxr * x_{s}(s=0)
% eq_PDE = diff(X,t)==lam*X+diff(X,s,2)+Bpv*x; % 	PDE: x_{t}  = lam * x + x_{ss} + Bpv * xo
% eq_BC = [subs(X,s,0)==0;subs(X,s,1)==0]; %  BCs: x(s=a) = 0,     x(s=b) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_ODE;eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
