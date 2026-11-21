% STABILITY_36 declares the 2D telegraph equation (below) on [0,1]x[0,1],
% displays it, converts it to a PIE and displays the PIE. Workspace: PDE
% ('pde_struct') and PIE ('pie_struct'). If a cell array params of
% assignment strings (e.g. {'lam=2;'}) exists in the workspace, its
% entries are evaluated after the default parameters.
%
%   Telegraph equation:
%   PDE: x1_{t} = x2
%        x2_{t} = -2*lam*x2 + C^2*(x1_{s1s1} + x1_{s2s2})
%   BCs: x1 = 0 and x2 = 0 on s1=0, s1=1, s2=0, and s2=1
%   C = 1; lam = 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_36
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
% Initial coding DB - 05/31/2026

% Clear the workspace
clc; clearvars -except params; clear stateNameGenerator
pvar s1 s2 t; % define independent variables
%% Specify parameters
C = 1;
lam = 1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable
x1 = pde_var(1,[s1;s2],[0,1;0,1]);
x2 = pde_var(1,[s1;s2],[0,1;0,1],[2;2]);

%% Define equations
eq_PDE = [diff(x1,'t')==x2;
          diff(x2,'t')==-2*lam*x2+C^2*(diff(x1,s1,2)+diff(x1,s2,2))];
eq_BC = [subs([x1;x2],s1,0)==0;
         subs([x1;x2],s1,1)==0;
         subs([x1;x2],s2,0)==0;
         subs([x1;x2],s2,1)==0];

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
