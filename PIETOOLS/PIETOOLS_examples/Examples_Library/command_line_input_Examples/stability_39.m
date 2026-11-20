% STABILITY_39 declares the linearized 2D isentropic compressible
% Navier-Stokes equations (below) on [0,1]x[0,1], with M = 0.1, lam = 1,
% nu = 1; displays them, converts them to a PIE and displays the PIE.
% Workspace: PDE ('pde_struct') and PIE ('pie_struct'). If a cell array
% params of assignment strings (e.g. {'M=0.2;'}) exists in the workspace,
% its entries are evaluated after the default parameters.
%
%   Linearized 2D isentropic compressible Navier-Stokes:
%   p_{t}  = -s2*p_{s1} - v1_{s1} - v2_{s2}
%   v1_{t} = -s2*v1_{s1} - v2 - (1/M^2)*p_{s1}
%            + nu*(v1_{s1s1}+v1_{s2s2}) + lam*(v1_{s1s1}+v2_{s2s1})
%   v2_{t} = -s2*v2_{s1} - (1/M^2)*p_{s2}
%            + nu*(v2_{s1s1}+v2_{s2s2}) + lam*(v1_{s1s2}+v2_{s2s2})
%   BCs: p(s1=0)=0, p(s2=0)=0, and v=0 on all edges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - stability_39
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
clc; clear; clear stateNameGenerator
pvar s1 s2 t; % define independent variables
%% Specify parameters
M = 0.1;
lam = 1;
nu = 1;
if exist('params','var')
    npars = length(params);
    %%% Specify potential parameters
    for j=1:npars
        eval(params{j});
    end
end
%% Define dependent variables and system variable
v1 = pde_var([s1;s2],[0,1;0,1]);
v2 = pde_var([s1;s2],[0,1;0,1]);
p = pde_var([s1;s2],[0,1;0,1]);

%% Define equations
eq_PDE = [diff(p,'t')==-s2*diff(p,s1)-diff(v1,s1)-diff(v2,s2);
          diff(v1,'t')==-s2*diff(v1,s1)-v2-(1/M^2)*diff(p,s1) ...
              +nu*(diff(v1,s1,2)+diff(v1,s2,2)) ...
              +lam*(diff(v1,s1,2)+diff(v2,[s1;s2]));
          diff(v2,'t')==-s2*diff(v2,s1)-(1/M^2)*diff(p,s2) ...
              +nu*(diff(v2,s1,2)+diff(v2,s2,2)) ...
              +lam*(diff(v1,[s1;s2])+diff(v2,s2,2))];
eq_BC = [subs([v1;v2],s1,0)==0;
         subs([v1;v2],s1,1)==0;
         subs([v1;v2],s2,0)==0;
         subs([v1;v2],s2,1)==0;
         subs(p,s1,0)==0;
         subs(p,s2,0)==0];

%% initialize pde system;
PDE = initialize([eq_PDE; eq_BC]);
%% display pde to verify and convert to pie
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
