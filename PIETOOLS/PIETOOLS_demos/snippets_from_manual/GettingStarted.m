
% % This script illustrates how PIETOOLS can be used to test stability of
% % a simple reaction-diffusion PDE:
% %
% % \dot{x}(t,s) = d^2/ds^2 x(t,s) + lam * x(t,s)              s in [0,1]
% %       x(t,0) = x(t,1) = 0
% %
% % with lam = 9. The script declares the PDE with pde_var, converts it to
% % a PIE, and runs the stability executive with the 'veryheavy' settings,
% % which prints the solver output. Workspace: pde, pie, the solved program
% % prog and the Lyapunov operator Pop. ECHO ON prints each line as it runs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - GettingStarted
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
% Initial coding DJ - 12/15/2022
% DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the stability test runs
%   lpiscript_sop; help text, license and authors added.
% MMP, 11/28/2026: _sop dropped from the names of the executives, LPI
%   and settings functions (2025 names where they exist); names only.

%% Clear the workspace; print each line as it runs
clc; clear;
echo on
%%%%%%%%%%%%%%%%%% Start Code Snippet %%%%%%%%%%%%%%%%%%
% Declare spatial and temporal variables s and t.
pvar s t; 
% % % Declare the PDE
% Declare state variable x.
x = pde_var('state',1,s,[0,1]);
% Declare the value of reaction parameter lam.
lam = 9;
% Declare the dynamics: d/dt x(t,s) = d^2/ds^2 x(t,s) + lam * x(t,s);
dynamics = (diff(x,t) == diff(x,s,2) + lam*x);
% Declare the boundary conditions: x(t,s=0) = x(t,s=1) = 0;
BC = [subs(x,s,0)==0; subs(x,s,1)==0];
% Add the dynamics and BCs to the PDE system
pde = initialize([dynamics;BC]);


% % % Convert the PDE to a PIE
pie = convert(pde,'pie');


% % % Run LPI stability test for PIE, using 'veryheavy' settings
%   [prog, Pop] = lpiscript(pie,  'stability','veryheavy');                 % DJ, MMP, 11/21/2026 (was)
[prog, Pop] = lpiscript(pie,  'stability','veryheavy');                     % DJ, MMP, 11/21/2026

%%%%%%%%%%%%%%%%%% End Code Snippet %%%%%%%%%%%%%%%%%%
echo off