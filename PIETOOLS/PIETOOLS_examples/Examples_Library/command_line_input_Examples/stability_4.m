% STABILITY_4 declares the 2x2 hyperbolic system (below) on [0,1] with
% r1 = 0.8, r2 = 1.1, sig1 = 2.3, sig2 = -3.5, qb = -0.7, pb = 0.5;
% displays it, converts it to a PIE and displays the PIE. Workspace: PDE
% ('pde_struct') and PIE ('pie_struct'). If a cell array params of
% assignment strings (e.g. {'qb=-0.5;'}) exists in the workspace, its
% entries are evaluated after the default parameters.
%
% 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s}           | Different parameters          Saba 2019 [3] 
%        x2_{t} = sig2*x1 + (1/r2)*x2_{s}           | may be invoked calling  
%   BCs: x1(s=0) = qb*x2(s=0)                       | examples 4.1, 4.2, 
%        x2(s=1) = pb*x1(s=1)                       | 4.3, 4.4, and 4.5.  
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable ) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_4
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
% Initial coding SS - 11/16/2022
% DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
%   workspace except a cell params of assignment strings, whose entries
%   are evaluated after the default parameters; the BCs read x1(0) =
%   qb*x2(0), x2(1) = pb*x1(1) (were x1(0) = qb*x2(1), x2(1) = pb*x1(0));
%   the results are named PDE and PIE (were pde and pie); help and
%   comments revised.

% Clear the workspace
%clc; clear; clear stateNameGenerator                                       % DB, MMP, 11/21/2026 (was)
clc; clearvars -except params; clear stateNameGenerator                     % DB, MMP, 11/21/2026
pvar s t; % define independent variables
% Define dependent variables and system variable
x = pde_var(2,s,[0,1]); % x = [x1; x2];
%% Specify the parameters
r1=.8;   r2=1.1;   sig1=2.3;   sig2=-3.5;    qb=-.7;   pb=.5; 
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==[0,sig1;sig2,0]*x+[-1/r1,0;0,1/r2]*diff(x,s); % 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s},x2_{t} = sig2*x1 + (1/r2)*x2_{s}
B1 = [[1,0]*subs(x,s,0);
      [0,1]*subs(x,s,1)];
B2 = [[0,1]*subs(x,s,0);
      [1,0]*subs(x,s,1)];
eq_BC = B1 == [qb,0; 0,pb] * B2;
%   BCs: x1(s=0) = qb*x2(s=0), x2(s=1) = pb*x1(s=1)                      
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% % 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s}           | Different parameters          Saba 2019 [3] 
% %        x2_{t} = sig2*x1 + (1/r2)*x2_{s}           | may be invoked calling  
% %   BCs: x1(s=0) = qb*x2(s=0)                       | examples 4.1, 4.2, 
% %        x2(s=1) = pb*x1(s=1)                       | 4.3, 4.4, and 4.5.
% x = state('pde',2); % x = [x1; x2];
% pde = sys();
% r1=.8;   r2=1.1;   sig1=2.3;   sig2=-3.5;    qb=-.7;   pb=.5; 
% %% Define equations
% eq_PDE = diff(x,t)==[0,sig1;sig2,0]*x+[-1/r1,0;0,1/r2]*diff(x,s); % 	PDE: x1_{t} = sig1*x2 - (1/r1)*x1_{s},x2_{t} = sig2*x1 + (1/r2)*x2_{s}
% eq_BC = [[1,0]*subs(x,s,0); [0,1]*subs(x,s,1)]==[0,qb;pb,0]*[[1,0]*subs(x,s,0);[0,1]*subs(x,s,1)];    
% %   BCs: x1(s=0) = qb*x2(s=0), x2(s=1) = pb*x1(s=1)                      
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
