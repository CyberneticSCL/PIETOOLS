% STABILITY_2 declares the 2x2 hyperbolic system x = [x_-; x_+] on [0,1],
%   x_t = Fm*x - Lm*x_s,
%   [x_-(1); x_+(0)] = [Gm1, Gm2; Gm3, Gm4]*[x_-(0); x_+(1)],
% with Fm, Lm and Gm1, ..., Gm4 as below; displays it, converts it to a
% PIE and displays the PIE. Workspace: PDE ('pde_struct') and PIE
% ('pie_struct'). If a cell array params of assignment strings (e.g.
% {'Gm1=0.5;'}) exists in the workspace, its entries are evaluated after
% the default parameters.
%
%	PDE: x_{t} = Fm*x - Lm*x_{s}                    | Different parameters          Lamare 2016 [1]
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]         | may be invoked calling
%        [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]           | examples 2.1, 2.2, 
%                                                                            | 2.3, and 2.4.    
%                                                                         (Exponentially PDE stable )
%                                                                         (Finite Energy PDE stable )   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_2
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
% MMP, 11/19/2026: Release comment pass (comment-only).

% Clear the workspace
clc; clear; clear stateNameGenerator
pvar s t; % define independent variables
%% Define dependent variables and system variable
%	PDE: x_{t} = Fm*x - Lm*x_{s}                    | Different parameters          Lamare 2016 [1] 
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]         | may be invoked calling
%        [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]         | examples 2.1, 2.2, 
%                                                   | 2.3, and 2.4.
x = pde_var(2,s,[0,1]); % x = [x_-; x_+];
%% Specify the parameters
Gm1=[.2]; Gm2=[-.3]; Gm3=[.6]; Gm4=[.1]; Lm=[-3 0;0 1]; Fm=[.2 -.3; .6 .1];
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define equations
eq_PDE = diff(x,t)==Fm*x-Lm*diff(x,s); % 	PDE: x_{t} = Fm*x - Lm*x_{s}
eq_BC = [[1,0]*subs(x,s,1); [0,1]*subs(x,s,0)]==[Gm1,Gm2;Gm3,Gm4]*[[1,0]*subs(x,s,0);[0,1]*subs(x,s,1)];    
%   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)], [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]
%% initialize pde system;
PDE = initialize([eq_PDE;eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
% %% Define dependent variables and system variable
% %	PDE: x_{t} = Fm*x - Lm*x_{s}                    | Different parameters          Lamare 2016 [1] 
% %   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)]         | may be invoked calling
% %        [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]         | examples 2.1, 2.2, 
% %                                                   | 2.3, and 2.4.
% x = state('pde',2); % x = [x_-; x_+];
% pde = sys();
% Gm1=[.2]; Gm2=[-.3]; Gm3=[.6]; Gm4=[.1]; Lm=[-3 0;0 1]; Fm=[.2 -.3; .6 .1];
% %% Define equations
% eq_PDE = diff(x,t)==Fm*x-Lm*diff(x,s); % 	PDE: x_{t} = Fm*x - Lm*x_{s}
% eq_BC = [[1,0]*subs(x,s,1); [0,1]*subs(x,s,0)]==[Gm1,Gm2;Gm3,Gm4]*[[1,0]*subs(x,s,0);[0,1]*subs(x,s,1)];    
% %   BCs: [x_-(s=1)] = [Gm1, Gm2] [x_-(s=0)], [x_+(s=0)] = [Gm3, Gm4] [x_+(s=1)]
% %% addequations to pde system; set control inputs/observed inputs, if any
% pde = addequation(pde,[eq_PDE;eq_BC]);
% %% display pde to verify and convert to pie
% display(pde);
% pie = convert(pde,'pie');
% display(pie);
