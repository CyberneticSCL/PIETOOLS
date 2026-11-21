% DEMO6_Hinf_optimal_control.m
%
% This document illustrates how an Hinfty optimal controller can be designed
% for a PDE using the PIE/LPI framework.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + s*w(t) + s*u(t),   s in [0,1];
% Outputs:       z(t)    = [int_{0}^{1} x(t,s) ds + w(t); u(t)];
% BCs:                0  = x(t,0) = d/ds x(t,1);
% (unstable for lam > pi^2/4)
% Letting v:=(d^2/ds^2)x, We derive an equivalent PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s) + [B2 u](t,s);
%               z(t) = [C1 v](t)  + [D11 w](t) + [D12 u](t);
% Using a state feedback control u = K*v we get the closed loop PIE
%   [T \dot{v}](t,s) = [(A + B2*K) v](t,s) + [B1 w](t,s);
%               z(t) = [(C1 + D12*K) v](t) + [D11 w](t)
% To compute an operator K that minimizes the L2 gain from disturbance w to
% the output z, we solve the LPI
%   min_{gam,P,Z}   gam
%   s.t.            P>=0,
%       [ -gam*I            D11      (C1*P+D12*Z)*T'            ]
%       [  D11'            -gam*I    B1'                        ] <= 0
%       [  T*(C1*P+D12*Z)'  B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
%
% Then, using K = Z*P^{-1}, the closed-loop system with u = K*v satisfies
%   ||z||_{L2}/||w||_{L2} <= gam
% We manually declare this LPI here; the commented call to lpiscript_sop
% after the closed-loop construction solves it with the predefined
% 'hinf-controller' executive instead.
% We simulate the open loop and closed loop response of the PDE for various
% initial conditions using PIESIM.
%
% The script displays the PDE (lam = 5), solves the LPI with P >= 1e-3*I,
% computes K from the solved P and Z (getController_sop), builds the
% closed-loop PIE PIE_CL, simulates the open and closed loop from
% x(0,s) = sin(pi*s/2) with w(t) = sin(pi*t)/t, and plots the state,
% output and control effort. Workspace: gam_val, Pval, Zval, Kval,
% PIE_CL, solution_OL, solution_CL and the sampled x_OL, z_OL, x_CL,
% z_CL, u_CL. ECHO ON prints each line as it runs.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO6
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
% MP, SS, DJ, 2022: Initial coding;
% DJ, 10/20/2024: Update to use new LPI programming functions;
% DJ, 11/19/2024: Simplify demo (remove lines of code where possible, get 
%                   rid of extra simulations);
% DJ, 12/15/2024: Use PIESIM_plotsolution to plot simulation results;
% DJ, 12/22/2024: Use piess;
% DB, 12/29/2024: Use pde_var objects instead of sys and state;
% YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: the LPI is declared
%   by lpiprogram(PIE.vars,PIE.dom), P = poscopvar(prog,T',4) +
%   1e-3*eye_copvar_sop(P), Z = lpivar_cdopvar(prog,B2',2) and
%   lpi_ineq_sop, in place of poslpivar, lpivar and lpi_ineq; lpigetsol_sop
%   reads the solution and getController_sop computes K; the closed loop is
%   a pie_struct; the predefined alternative is lpiscript_sop; PIESIM takes
%   uinput.ic as an array, no ndiff, and returns solution.timedep.primary
%   and regulated as cells; help text revised.

% Clear the workspace and figures; print each line as it runs
clear; clc; close all; clear stateNameGenerator
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar t s
% Declare state, input, and output variables
x = pde_var('state',1,s,[0,1]);   
z = pde_var('output',2);
w = pde_var('input',1);          u = pde_var('control',1);
% Declare the sytem equations
lam = 5;
PDE = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
       z == [int(x,s,[0,1])+w; u];
       subs(x,s,0)==0;
       subs(diff(x,s),s,1)==0];

% Display the PDE
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE,'pie');
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(PIE.vars,PIE.dom);                                        % YP, DJ, MMP, 11/21/2026

% % Declare decision variables:
% %   gam \in \R,     P:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable P>=0
[prog,P] = poscopvar(prog,T',4);        % P on the input spaces of T        % YP, DJ, MMP, 11/21/2026
% Enforce strict positivity P >= 1e-3
P = P + 1e-3*eye_copvar_sop(P);                                             % YP, DJ, MMP, 11/21/2026
% Indefinite operator variable Z: L2-->\R
[prog,Z] = lpivar_cdopvar(prog,B2',2);  % Z with the spaces of B2'          % YP, DJ, MMP, 11/21/2026

% % Set inequality constraints:
% %   Q <= 0
Q = [-gam*eye(2),       D11,    (C1*P+D12*Z)*(T');
     D11',              -gam,   B1';
     T*(C1*P+D12*Z)',   B1,     (A*P+B2*Z)*(T')+T*(A*P+B2*Z)'];
prog = lpi_ineq_sop(prog,-Q);                                               % YP, DJ, MMP, 11/21/2026

% % Set objective function:
% %   min gam
prog = lpisetobj(prog, gam);

% % Solve and retrieve the solution
prog_sol = lpisolve(prog);
% Extract solved value of decision variables
gam_val = double(lpigetsol_sop(prog_sol,gam));                              % YP, DJ, MMP, 11/21/2026
Pval = lpigetsol_sop(prog_sol,P);                                           % YP, DJ, MMP, 11/21/2026
Zval = lpigetsol_sop(prog_sol,Z);                                           % YP, DJ, MMP, 11/21/2026
% Build the optimal control operator K.
[Kval,~,~] = getController_sop(Pval,Zval);                                  % YP, DJ, MMP, 11/21/2026

% % Construct the closed-loop PIE system
%   d/dt T*v(t,s) = (A + B2*K)*v(t,s) + B1*w(t,s);
%            z(t) = (C1 + D12*K)*v(t) + D11*w(t)
PIE_CL = pie_struct();                                                      % YP, DJ, MMP, 11/21/2026
PIE_CL.T = T;     PIE_CL.A = A+B2*Kval;     PIE_CL.Bw = B1;                 % YP, DJ, MMP, 11/21/2026
PIE_CL.Cz = C1+D12*Kval;     PIE_CL.Dzw = D11;                              % YP, DJ, MMP, 11/21/2026
PIE_CL = initialize(PIE_CL);                                                % YP, DJ, MMP, 11/21/2026

% % Alternatively, uncomment to use pre-defined functions
% [prog, Kval, gam_val] = lpiscript_sop(PIE,'hinf-controller','heavy');     % YP, DJ, MMP, 11/21/2026
% PIE_CL = closedLoopPIE(PIE,Kval);


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput.ic = sin(sx*pi/2);
uinput.w = sin(pi*st)./(st+eps); 

% % Set options for discretization and simulation
opts.plot = 'yes';   % plot the solution
opts.N = 16;        % expand using 16 Chebyshev polynomials
opts.tf = 2;        % simulate up to t = 2
opts.dt = 1e-2;     % use time step of 10^-2

% % Perform the actual simulation
% Simulate uncontrolled PIE and extract solution
[solution_OL,grid] = PIESIM(PIE,opts,uinput);
tval = solution_OL.timedep.dtime;
x_OL = reshape(solution_OL.timedep.primary{2}(:,1,:),opts.N+1,[]);
z_OL = solution_OL.timedep.regulated{1}(1,:);
% Simulate controlled PIE and extract solution
[solution_CL,~] = PIESIM(PIE_CL,opts,uinput);
x_CL = reshape(solution_CL.timedep.primary{2}(:,1,:),opts.N+1,[]);
z_CL = solution_CL.timedep.regulated{1}(1,:);
u_CL = solution_CL.timedep.regulated{1}(2,:);
w = double(subs(uinput.w,st,tval));


echo off


% % Change titles of plots
fig1 = figure(1);       ax1 = fig1.CurrentAxes;
title(ax1,'Open-Loop Regulated Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',15);
fig2 = figure(2);
fig2.Children.Title.String = ['Open-Loop Final PDE State $\mathbf{x}(t=',num2str(opts.tf),',s)$'];
fig3 = figure(3);
fig3.Children(2).Title.String = 'Open-Loop PDE State Evolution $\mathbf{x}(t,s)$';
fig4 = figure(4);       ax4 = fig4.CurrentAxes;
title(ax4,'Closed-Loop Regulated Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',15);
fig5 = figure(5);
fig5.Children.Title.String = ['Closed-Loop Final PDE State $\mathbf{x}(t=',num2str(opts.tf),',s)$'];
fig6 = figure(6);
fig6.Children(2).Title.String = 'Closed-Loop PDE State Evolution $\mathbf{x}(t,s)$';
