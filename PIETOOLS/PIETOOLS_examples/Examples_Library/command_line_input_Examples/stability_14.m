% STABILITY_14 declares the Timoshenko beam in the hyperbolic states
% x1, ..., x4 on [0,1]:
%   x1_t = (1/aa)*x2_s,   x2_t = k*aa*g*(x1_s - x3),
%   x3_t = (1/II)*(x2 + x4_s) - c*x3,   x4_t = II*x3_s,
%   x1(0) = x3(0) = 0,   x2(1) = x4(1) = 0,
% with k = aa = II = g = 1 and damping c = 0; displays it, converts it to
% a PIE and displays the PIE. Workspace: PDE ('pde_struct') and PIE
% ('pie_struct'). If a cell array params of assignment strings (e.g.
% {'c=0.1;'}) exists in the workspace, its entries are evaluated after the
% default parameters.
%
% Timoshenko Beam
%   PDE: r*aa * w_{tt} = k*aa*g * (-phi_{s} + w_{ss})                               Peet 2019 [8] (Example 8.1.0.2)
%        r*II * phi_{tt} = E*II * phi_{ss}  + k*aa*g * (w_{s} - phi)-c
%        phi_{t}  
%   BCs: phi(s=0) = 0,      phi_{s}(s=1) = 0    
%        w(s=0) = 0,        w_{s}(s=1) - phi(s=1) = 0           
% 
%   Use states:                                     | k = 1                 Hyperbolic implementation | (PDE stable without damping c=0)
%        x1 = w_{t},   x2 = k*aa*g * (w_{s}-phi),   | aa = 1 
%        x3 = phi_{t}-c phi_{t}, x4 = E*II * phi_{s}.         | II = 1   
%       =>                                          | g = 1
%   PDE: x1_{t} = (1/r/aa) * x2_{s}                 | E = 1   (PDE stable and Finite-energy PDE stable with damping c >0)
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3     | r = 1
%        x3_{t} = (1/r/II)*x2 + (1/r/II)*x4_{s}     | 
%        x4_{t} = E*II * x3_{s}                     |
%   BCs: x1(0) = 0,         x2(1) = 0               |
%        x3(0) = 0,         x4(1) = 0               |
%                                                   |

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_14
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
% DB, MMP, 11/21/2026: Changes since PIETOOLS 2025: the script clears the
%   workspace except a cell params of assignment strings, whose entries
%   are evaluated after the default parameters; new parameter c (default
%   0), the damping term -c*x3 in the dynamics of x3; the results are
%   named PDE and PIE (were pde and pie); help and comments revised.

% Clear the workspace
%clc; clear; clear stateNameGenerator                                       % DB, MMP, 11/21/2026 (was)
clc; clearvars -except params; clear stateNameGenerator                     % DB, MMP, 11/21/2026
pvar s t; % define independent variables
%% Specify the parameters 
k = 1;    aa = 1;    II = 1;    g = 1;    c=0;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable    
pde_var x1 x2 x3 x4
x1.vars = s;    x2.vars = s;    x3.vars = s;    x4.vars = s;
%% Define equations
%   PDE: x1_{t} = (1/aa) * x2_{s}                
%        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
%        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}  -c x_3   
%        x4_{t} = II * x3_{s}
eq_PDE = [diff(x1,t)==(1/aa)*diff(x2,s); 
          diff(x2,t)==k*aa*g*diff(x1,s)-k*aa*g*x3;
          diff(x3,t)==(1/II)*x2+(1/II)*diff(x4,s)-c*x3;
          diff(x4,t)==II*diff(x3,s)];
eq_BC = [subs(x1,s,0)==0;subs(x3,s,0)==0;
         subs(x2,s,1)==0;subs(x4,s,1)==0]; %   BCs:BCs: x1(0) = 0, x2(1) = 0, x3(0) = 0, x4(1) = 0 
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);



% %% Define dependent variables and system variable
% %   PDE: x1_{t} = (1/aa) * x2_{s}                
% %        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
% %        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}     
% %        x4_{t} = II * x3_{s}                   
% %   BCs: x1(0) = 0,         x2(1) = 0           
% %        x3(0) = 0,         x4(1) = 0           
% x1=state('pde');x2=state('pde');x3=state('pde');x4=state('pde'); 
% pde = sys(); 
% k = 1;    aa = 1;    II = 1;    g = 1;    
% %% Define equations
% %   PDE: x1_{t} = (1/aa) * x2_{s}                
% %        x2_{t} = k*aa*g * x1_{s} - k*aa*g * x3   
% %        x3_{t} = (1/II)*x2 + (1/II)*x4_{s}     
% %        x4_{t} = II * x3_{s}
% eq_PDE = [diff(x1,t)==(1/aa)*diff(x2,s); 
%           diff(x2,t)==k*aa*g*diff(x1,s)-k*aa*g*x3;
%           diff(x3,t)==(1/II)*x2+(1/II)*diff(x4,s);
%           diff(x4,t)==II*diff(x3,s)];
% eq_BC = [subs(x1,s,0)==0;subs(x3,s,0)==0;
%          subs(x2,s,1)==0;subs(x4,s,1)==0]; %   BCs:BCs: x1(0) = 0, x2(1) = 0, x3(0) = 0, x4(1) = 0 
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
