% DEMO6_Hinf_optimal_control.m
% See Chapter 11.6 of the manual for a description.
%
% This document illustrates how an Hinfty optimal controller can be designed
% for a PDE using the PIE/LPI framework.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t) + u(t),   s in [0,1];
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
%       [  T*(C1*P+D12*Z)   B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
%
% Then, using K = Z*P^{-1}, the closed-loop system with u = K*v satisfies
%   ||z||_{L2}/||w||_{L2} <= gam
% We manually declare this LPI here, but it can also be solved using the
% "PIETOOLS_Hinf_control" executive file.
% We simulate the open loop and closed loop response of the PDE for various
% initial conditions using PIESIM.
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
% MP, SS, DJ, 2022: Initial coding;
% DJ, 10/20/2024: Update to use new LPI programming functions;
% DJ, 11/19/2024: Simplify demo (remove lines of code where possible, get 
%                   rid of extra simulations);
% DJ, 12/15/2024: Use PIESIM_plotsolution to plot simulation results;
% DJ, 12/22/2024: Use piess;

clc; clear; close all;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar t s
% Declare state, input, and output variables
x = state('pde');       w = state('in');
z = state('out', 2);    u = state('in');
% Declare the sytem equations
lam = 5;
PDE = sys();
eqs = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
       z == [int(x,s,[0,1])+w; u];
       subs(x,s,0)==0;
       subs(diff(x,s),s,1)==0];
PDE = addequation(PDE,eqs);
% Set u as constrolled input
PDE = setControl(PDE,u);
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE,'pie');
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    B2 = PIE.B2;
B1 = PIE.B1;    D11 = PIE.D11;  D12 = PIE.D12;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(s,[0,1]);

% % Declare decision variables:
% %   gam \in \R,     P:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable P>=0
[prog,P] = poslpivar(prog,[0,0;1,1],4);
% Enforce strict positivity P >= 1e-3
P = P + 1e-3;
% Indefinite operator variable Z
[prog,Z] = lpivar(prog,[1,0;0,1],2);

% % Set inequality constraints:
% %   Q <= 0
Q = [-gam*eye(2),       D11,    (C1*P+D12*Z)*(T');
     D11',              -gam,   B1';
     T*(C1*P+D12*Z)',   B1,     (A*P+B2*Z)*(T')+T*(A*P+B2*Z)'];
prog = lpi_ineq(prog,-Q);

% % Set objective function:
% %   min gam
prog = lpisetobj(prog, gam);

% % Solve and retrieve the solution
opts.solver = 'sedumi';         % Use SeDuMi to solve the SDP
opts.simplify = true;           % Simplify the SDP before solving
prog_sol = lpisolve(prog,opts);
% Extract solved value of decision variables
gam_val = lpigetsol(prog_sol,gam);
Pval = lpigetsol(prog_sol,P);
Zval = lpigetsol(prog_sol,Z);
% Build the optimal control operator K.
Kval = getController(Pval,Zval,1e-3);

% % Construct the closed-loop PIE system
%   d/dt T*v(t,s) = (A + B2*K)*v(t,s) + B1*w(t,s);
%            z(t) = (C1 + D12*K)*v(t) + D11*w(t)
PIE_CL = piess(T,A+B2*Kval,B1,C1+D12*Kval,D11);

% % Alternatively, uncomment to use pre-defined functions
% [prog, Kval, gam_val] = lpiscript(PIE,'hinf-controller','heavy');   
% PIE_CL = closedLoopPIE(PIE,Kval);


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput.ic.PDE = sin(sx*pi/2);
uinput.w = sin(pi*st)./(st+eps); 

% % Set options for discretization and simulation
opts.plot = 'no';   % don't plot final solution
opts.N = 16;        % expand using 16 Chebyshev polynomials
opts.tf = 2;        % simulate up to t = 2
opts.dt = 1e-2;     % use time step of 10^-2
ndiff = [0,0,1];    % PDE state involves 1 second order differentiable state variable

% % Perform the actual simulation
% Simulate uncontrolled PIE and extract solution
[solution_OL,grid] = PIESIM(PIE,opts,uinput,ndiff);
tval = solution_OL.timedep.dtime;
x_OL = reshape(solution_OL.timedep.pde(:,1,:),opts.N+1,[]);
z_OL = solution_OL.timedep.regulated(1,:);
% Simulate controlled PIE and extract solution
[solution_CL,~] = PIESIM(PIE_CL,opts,uinput,ndiff);
x_CL = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
z_CL = solution_CL.timedep.regulated(1,:);
u_CL = solution_CL.timedep.regulated(2,:);
w = double(subs(uinput.w,st,tval));


echo off


% % Plot simulated states and regulated outputs against time.
figs_OL = PIESIM_plotsolution(solution_OL,grid,'title','Open-Loop');
figs_CL = PIESIM_plotsolution(solution_CL,grid,'title','Closed-Loop');
% Change titles of plots
fig2 = figs_OL{2};  ax2 = fig2.CurrentAxes;
title(ax2,'Open-Loop Regulated Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',15);
fig4 = figs_CL{2};  ax4 = fig4.CurrentAxes;
title(ax4,'Closed-Loop Regulated Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',15);