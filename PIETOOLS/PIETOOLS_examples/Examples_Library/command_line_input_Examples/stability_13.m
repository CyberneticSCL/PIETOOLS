% STABILITY_13 declares the Euler-Bernoulli beam u_tt = -c*u_ssss on [0,1]
% (c = 0.1) in the states x1 = u_t, x2 = u_ss (below); displays it,
% converts it to a PIE and displays the PIE. Workspace: PDE ('pde_struct')
% and PIE ('pie_struct'). If a cell array params of assignment strings
% (e.g. {'c=1;'}) exists in the workspace, its entries are evaluated after
% the default parameters.
%
%% Euler-Bernoulli Equation
% Parabolic/Diffusive implementation
% % This case uses PDE states x1 = u_{t}, x2 = u_{ss}, so that the system becomes
% % PDE        x1_{t} = -c * x2_{ss}
% %            x2_{t} = x1_{ss}                          | (PDE stable)
% % with BCs   x1(0) = 0, x2(1) = 0, x1_{s}(0) = 0, x2_{s}(1) = 0   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_13
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
%   params before it was read; eq_PDE ends in ';' (it was printed).

% Clear the workspace
%clc; clear; clear stateNameGenerator                                       % MMP, 11/20/2026 (was)
clc; clearvars -except params; clear stateNameGenerator                     % MMP, 11/20/2026
pvar s t; % define independent variables
% Define dependent variables and system variable
x1=pde_var(s,[0,1]); x2=pde_var(s,[0,1]); 
%% Specify the parameters
c = 0.1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = [diff(x1,t,1)==-c*diff(x2,s,2); % PDE        x1_{t} = -c * x2_{ss}
%                  diff(x2,t,1)==diff(x1,s,2);]%           x2_{t} = x1_{ss} % MMP, 11/20/2026 (was)
                   diff(x2,t,1)==diff(x1,s,2);];%           x2_{t} = x1_{ss} % MMP, 11/20/2026
eq_BC = [subs(x1,s,0)==0;subs(diff(x1,s,1),s,0)==0;% with BCs   x1(0) = 0, x1_{s}(0) = 0, 
         subs(x2,s,1)==0;subs(diff(x2,s,1),s,1)==0];  % x2(1) = 0, x2_{s}(1) = 0
%% initialize pde system; PIE states (u_{tss},u_{ssss})
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: u_{tt} = -c*u_{ssss}                       | c = 0.1
% %   BCs: u(s=0) = 0,        u_{ss}(s=1) = 0         |                               Peet 2019 [8] (Example 8.1.0.1)
% %        u_{s}(s=0) = 0,    u_{sss}(s=1) = 0        |  
% x=state('pde'); 
% pde = sys(); 
% c = 0.1;
% %% Define equations
% eq_PDE = diff(x,t,2)==-c*diff(x,s,4); % 	PDE: u_{tt} = -c*u_{ssss}
% eq_BC = [subs(x,s,0)==0;subs(diff(x,s),s,0)==0;
%          subs(diff(x,s,2),s,1)==0;subs(diff(x,s,3),s,1)==0]; %   BCs: u(s=0) = 0, u_{ss}(s=1) = 0, u_{s}(s=0) = 0, u_{sss}(s=1) = 0
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
