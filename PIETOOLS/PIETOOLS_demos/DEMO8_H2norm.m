% DEMO8_H2Norm.m
% See the manual for full description.
%
% This document illustrates, with a simple example, how PIETOOLS can be used to
% analyse and design optimal controllers and observers based on the H2 norm
%  The example is a pure convection equation in 1D: 
%  PDE :        x_{t} = x_{s} + (s-s^2)w(t)
%  With BC     x(s=1) = 0
%  And output  z(t) = int(x(t,s),s,0,1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO8
%
% Copyright (C)2024  PIETOOLS Team
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
% DB, 08/16/2024: Initial coding;
% DJ, 12/26/2024: Update to use new LPI programming functions;
% DJ, 12/15/2024: Match structure to that of other demos;
% DB, 12/27/2024: Replace sys() objects by pde_var()

clc; clear; close all;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = pde_var('state',1,s,[0,1]);
w = pde_var('in',1);
z = pde_var('out',1);
% Declare the sytem equations
pde = [diff(x,t,1)==diff(x,s,1)+(s-s^2)*w;    % dynamics
                z==int(x,s,[0,1]);                     % output equation
                subs(x,s,1)==0];                            % boundary condition
pde=initialize(pde);
display_PDE(pde);
% % Convert PDE to PIE
PIE = convert(pde,'pie');
T = PIE.T;      A = PIE.A;
B1 = PIE.B1;    C1 = PIE.C1;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(s,[0,1]); 

% % Declare decision variables:
% %   gam \in \R,     W:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable W>=0 with default polynomial degrees up to 3.
[prog,W] = poslpivar(prog,[0;1]);

% % Set inequality constraints:
% %   A W T* + T W A* + B1 B1* <= 0
% %   gam >= trace(C1 W C1*)
% Operator inequality Dop<=0
Dop =  A*W*T'+T*W*A'+B1*B1';
prog = lpi_ineq(prog,-Dop);
% Scalar inequality gam >= trace(C1 W C1*)
Aux = C1*W*C1';
traceVal = trace(Aux.P);
prog = lpi_ineq(prog, gam-traceVal);

% % Set objective function:
% %   min gam
prog = lpisetobj(prog, gam);

% % Solve and retrieve the solution
opts.solver = 'sedumi';         % Use SeDuMi to solve the SDP
prog_sol = lpisolve(prog,opts);
% Extract solved value of decision variables
gamd = sqrt(double(lpigetsol(prog_sol,gam)));
Wc = lpigetsol(prog_sol,W);


echo off

fprintf([' If succesful, ',num2str(gamd),' is an upper bound on the H2 norm of the specified system.\n']);