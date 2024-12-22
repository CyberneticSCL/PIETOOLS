% DEMO5_Hinf_optimal_estimator.m
% See Chapter 11.5 of the manual for a description.
%
% This document illustrates how an Hinfty optimal estimator can be designed
% for a PIE, and how the estimated state can be simulated.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t),   s in [0,1];
% Outputs:          z(t) = int_{0}^{1} x(t,s) ds + w(t);
%                   y(t) = x(t,1);
% BCs:                0  = x(t,0) = x_{s}(t,1);
% (unstable for lam > pi^2/4 = 2.4674)
% Letting v:=(d^2/ds^2)x, We derive an equivalent PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s);
%               z(t) = [C1 v](t)  + [D11 w](t);
%               y(t) = [C2 v](t)  + [D21 w](t);
% So that x = T*v. We design an estimator of the form:
%   [T \dot{vhat}](t,s) = [A vhat](t,s) + [L(yhat-y)](t,s);
%               zhat(t) = [C1 vhat](t)
%               yhat(t) = [C2 vhat](t)
% Then, the errors e=(vhat-v) and ztilde = (zhat-z) satisfy
%   [T \dot{e}](t,s) = [(A+L*C2) e](t,s) - [(B1+L*D21) w](t,s);
%          ztilde(t) = [C1 e](t)         - [D11 w](t);
% To compute an operator L that minimizes the L2 gain from
% disturbances w to error ztilde in the output, we solve the LPI
%   min_{gam,P,Z}   gam
%   s.t.            P>=0,
%       [-gam*I,           -D11',  -(P*B1+Z*D21)'*T           ]=: Q <=0
%       [-D11,             -gam*I, C1                         ]
%       [-T'*(P*B1+Z*D21), C1',    (P*A+Z*C2)'*T+T'*(P*A+Z*C2)]
%
% Then, using L = P^{-1}*Z, the L2 gain satisfies 
%   ||ztilde||_{L2}/||w||_{L2} <= gam
% We manually declare this LPI here, but it can also be solved using the
% "PIETOOLS_Hinf_estimator" executive file.
% We simulate the PDE state x and estimated state xhat using PIESIM.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO5
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
% DJ, 11/19/2024: Simplify demo (remove lines of code where possible);
% DJ, 12/15/2024: Use PIESIM_plotsolution to plot simulation results;
% DJ, 12/22/2024: Use piess and pielft;

clc; clear; close all;
echo on

% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar t s
% Declare state, input, and output variables
x = state('pde');   w = state('in');
y = state('out');   z = state('out');
% Declare the sytem equations
lam = 4;
PDE = sys();
eqs = [diff(x,t) == diff(x,s,2) + lam*x + w;    % PDE
       z == int(x,s,[0,1]) + w;                 % regulated output
       y == subs(x,s,1);                        % observed output
       subs(x,s,0) == 0;                        % first boundary condition
       subs(diff(x,s),s,1) == 0];               % second boundary condition
PDE = addequation(PDE,eqs);                     % Add the equations to the PDE structure
% Set y as an observed output
PDE = setObserve(PDE,y);                        
display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE);
T = PIE.T;
A = PIE.A;      C1 = PIE.C1;    C2 = PIE.C2;
B1 = PIE.B1;    D11 = PIE.D11;  D21 = PIE.D21;


% =============================================
% === Declare the LPI

% % Initialize LPI program
prog = lpiprogram(s,[0,1]);

% % Declare decision variables:
% %   gam \in \R,     P:L2-->L2,    Z:\R-->L2
% Scalar decision variable
[prog,gam] = lpidecvar(prog,'gam');
% Positive operator variable P>=0
opts.sep = 1;                   % set P.R.R1=P.R.R2
[prog,P] = poslpivar(prog,T.dim,4,opts);
% Indefinite operator variable Z: \R-->L2
[prog,Z] = lpivar(prog,[0,1;1,0],4);

% % Set inequality constraints:
% %   Q <= 0
Q = [-gam,              -D11',        -(P*B1+Z*D21)'*T;
     -D11,              -gam,         C1;
     -T'*(P*B1+Z*D21),  C1',          (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];
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
% Build optimal observer operator L
Lval = getObserver(Pval,Zval);

% % Construct the closed-loop PIE system
% Declare Luenberger estimator with input y
%   d/dt T*vhat(t,s) = (A+L*C2)*vhat(t,s) - Ly(t,s);
%            zhat(t) = C1*vhat(t);
% and take LFT with PIE
PIE_est = piess(T,A+Lval*C2,-Lval,C1(1,:));
PIE_CL = pielft(PIE,PIE_est);



% % Alternatively, uncomment to use pre-defined functions
% [prog, Lval, gam_val] = lpiscript(PIE,'hinf-observer','heavy');   
% PIE_CL = closedLoopPIE(PIE,Lval,'observer');


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput.ic.PDE = [-10*sx;    % actual initial PIE state value
                 0];        % estimated initial PIE state value
uinput.w = 2*sin(pi*st);

% % Set options for discretization and simulation
opts.plot = 'no';   % don't plot final solution
opts.N = 8;         % expand using 8 Chebyshev polynomials
opts.tf = 2;        % simulate up to t = 2
opts.dt = 1e-3;     % use time step of 10^-3
ndiff = [0,0,2];    % PDE state involves 2 second order differentiable state variables

% % Simulate solution to the PIE with estimator.
[solution,grid] = PIESIM(PIE_CL,opts,uinput,ndiff);
% % Extract actual and estimated state and output at each time step.
tval = solution.timedep.dtime;
x_act = reshape(solution.timedep.pde(:,1,:),opts.N+1,[]);
x_est = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
z_act = solution.timedep.regulated(1,:);
z_est = solution.timedep.regulated(2,:);


echo off


% % Plot simulated states and regulated outputs against time.
figs = PIESIM_plotsolution(solution,grid,'tlog');
% Change titles of plots
fig1 = figs{1};
fig1.Children(5).String = 'True ($x_1(t,s)$) and estimated ($x_2(t,s)$) PDE state evolution';
fig2 = figs{2};     ax = fig2.CurrentAxes;
title(ax,'True ($z_1(t)$) and estimated ($z_2(t)$) regulated output evolution','Interpreter','latex','FontSize',15);