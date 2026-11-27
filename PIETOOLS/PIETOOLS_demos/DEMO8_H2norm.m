% DEMO8_H2norm.m
%
% This document illustrates, with a simple example, how PIETOOLS can be used to
% compute an upper bound on the H2 norm of a PDE.
%  The example is a pure convection equation in 1D: 
%  PDE :        x_{t} = x_{s} + (s-s^2)w(t)
%  With BC     x(s=1) = 0
%  And output  z(t) = int(x(t,s),s,0,1)
%
% The script converts the PDE to a PIE with operators T, A, B1, C1 and
% solves the LPI
%   min gam   s.t.   W >= 0,   A*W*T' + T*W*A' + B1*B1' <= 0,
%                    gam >= trace(C1*W*C1'),
% so that sqrt(gam) bounds the H2 norm from w to z. It displays the PDE and
% prints the bound gamd = sqrt(gam). Workspace: the solved program
% prog_sol, gamd and the solved operator Wc. ECHO ON prints each line as it
% runs.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DB, 08/16/2024: Initial coding;
% DJ, 12/26/2024: Update to use new LPI programming functions;
% DJ, 12/15/2024: Match structure to that of other demos;
% DB, 12/27/2024: Replace sys() objects by pde_var();
% DJ, 05/18/2025: Declare no particular SDP solver by default;
% DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI is declared by
%   lpiprogram(PIE.vars,PIE.dom), W = poscopvar(prog,T',3) and lpi_ineq_sop,
%   with the trace of C1*W*C1' by trace_rn_sop, in place of poslpivar,
%   lpi_ineq and trace(Aux.P); lpigetsol_sop reads the solution; help text
%   revised.
% MMP, 11/27/2026: copvar, cdopvar renamed opvar, dopvar, with the
%   functions of the rename table; names only.


% Clear the workspace and figures; print each line as it runs
clc; clear; close all; clear stateNameGenerator;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar s t
% Declare state, input, and output variables
x = pde_var('state',s,[0,1]);
w = pde_var('input',1);
z = pde_var('output',1);
% Declare the sytem equations
pde = [diff(x,t,1)==diff(x,s,1)+(s-s^2)*w;    % dynamics
                 z==int(x,s,[0,1]);           % output equation
                subs(x,s,1)==0];              % boundary condition
% Initialize and display the PDE
pde=initialize(pde);
display_PDE(pde);
% % Convert PDE to PIE
PIE = convert(pde,'pie');
T = PIE.T;      A = PIE.A;
B1 = PIE.B1;    C1 = PIE.C1;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(PIE.vars,PIE.dom);                                        % DJ, MMP, 11/21/2026

% % Declare decision variables:
% %   gam \in \R,     W:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable W>=0 with default polynomial degrees up to 3.
[prog,W] = poslpivar(prog,T',3);        % W on the input spaces of T        % DJ, MMP, 11/21/2026

% % Set inequality constraints:
% %   A W T* + T W A* + B1 B1* <= 0
% %   gam >= trace(C1 W C1*)
% Operator inequality Dop<=0
Dop =  A*W*T'+T*W*A'+B1*B1';
prog = lpi_ineq_sop(prog,-Dop);                                             % DJ, MMP, 11/21/2026
% Scalar inequality gam >= trace(C1 W C1*)
Aux = C1*W*C1';
traceVal = trace_rn_sop(Aux);                                               % DJ, MMP, 11/21/2026
prog = lpi_ineq_sop(prog, gam-traceVal);                                    % DJ, MMP, 11/21/2026

% % Set objective function:
% %   min gam
prog = lpisetobj(prog, gam);

% % Solve and retrieve the solution
prog_sol = lpisolve(prog);
% Extract solved value of decision variables
gamd = sqrt(double(lpigetsol_sop(prog_sol,gam)));                           % DJ, MMP, 11/21/2026
Wc = lpigetsol_sop(prog_sol,W);                                             % DJ, MMP, 11/21/2026


echo off

% Print the bound on the H2 norm
fprintf([' If succesful, ',num2str(gamd),' is an upper bound on the H2 norm of the specified system.\n']);