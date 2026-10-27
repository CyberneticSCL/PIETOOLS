% DEMO7_observer_based_control.m
% See Chapter 11.7 of the manual for a description.
%
% This document illustrates construction of an observer-based state feedback 
% controller. We construct the observer and controller separately
% and couple them to find a closed loop system which is then simulated for 
% different IC and disturbances.
% Specifically, we consider the following system:
% PDE:      \dot{x}(t,s) = (d^2/ds^2) x(t,s) + lam x(t,s) + w(t) + u(t),   s in [0,1];
% Outputs:       z(t)    = [int_{0}^{1} x(t,s) ds; u(t)];
% BCs:                0  = x(t,0) = d/ds x(t,1);
% (unstable for lam > pi^2/4)
% First we convert the above PDE to a PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s) + [B1 w](t,s) + [B2 u](t,s);
%               z(t) = [C1 v](t)  + [D11 w](t) + [D12 u](t);
%               y(t) = [C2 v](t)  + [D21 w](t) + [D22 u](t)
% We find the observer gain L by solving the LPI
%   min_{gam,P,Z} gam
%   s.t.          P>=0,
%       [-gam*I,           -D11',  -(P*B1+Z*D21)'*T           ] <=0
%       [-D11,             -gam*I, C1                         ]
%       [-T'*(P*B1+Z*D21), C1',    (P*A+Z*C2)'*T+T'*(P*A+Z*C2)]
% where L = P^{-1}*Z.
% Likewise we find a state feedback control u = Kv by solving the LPI
%   min_{gam,P,Z} gam
%   s.t.          P>=0,
%       [ -gam*I            D11      (C1*P+D12*Z)*T'            ]
%       [  D11'            -gam*I    B1'                        ] <= 0
%       [  T*(C1*P+D12*Z)   B1       (A*P+B2*Z)*T'+T*(A*P+B2*Z)']
% where K = Z*P^{-1}. Then, constructing an estimator with gain L, and
% using feedback control u = K*hat{v}, the closed-loop PIE takes the form
%   [T, 0] [\dot{v}(t)   ] = [A,     B2*K  ] [v(t)   ] + [B1   ] w(t)
%   [0, T] [\dot{vhat}(t)]   [-L*C2, A+L*C2] [vhat(t)] + [L*D21]
%
%                [z(t)   ] = [C1, D12*K ] [v   ] + [D11] w(t)
%                [zhat(t)]   [0,  C1    ] [vhat]   [0  ]
% We simulate the open- and closed-loop system response for various
% initial conditions using PIESIM.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO7
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
% DJ, 12/22/2024: Use piess and pielft;
% DB, 12/29/2024: Use pde_var objects instead of sys and state

clear; clc; close all; clear stateNameGenerator
echo on


% =============================================
% === Declare the system of interest

% % Declare system as PDE
% Declare independent variables (time and space)
pvar t s
% Declare state, input, and output variables
x = pde_var('state',1,s,[0,1]);   
w = pde_var('input',1);      u = pde_var('control',1);
z = pde_var('output',2);     y = pde_var('sense',1);
% Declare the sytem equations
lam = 5;
PDE = [diff(x,t) == diff(x,s,2) + lam*x + s*w + s*u;
       z == [int(x,s,[0,1]); u];
       y == subs(x,s,1);
       subs(x,s,0)==0;
       subs(diff(x,s),s,1)==0];

display_PDE(PDE);

% % Convert PDE to PIE
PIE = convert(PDE,'pie');
T = PIE.T;      
A = PIE.A;      B1 = PIE.B1;    B2 = PIE.B2;
C1 = PIE.C1;    D11 = PIE.D11;  D12 = PIE.D12;
C2 = PIE.C2;    D21 = PIE.D21;  D22 = PIE.D22;


% =============================================
% === Declare the LPI

% % Use the predefined Hinf control and estimator functions.
settings = lpisettings('heavy');
[prog_k, Kval, gam_co_val] = lpiscript(PIE,'hinf-controller',settings);
[prog_l, Lval, gam_ob_val] = lpiscript(PIE,'hinf-observer',settings);

% % Construct the closed-loop PIE system
% Declare Luenberger estimator with input y and output u
%   d/dt T*vhat(t,s) = (A+L*C2)*vhat(t,s) - Ly(t,s);
%            zhat(t) = C1*vhat(t);
%               u(t) = K*vhat(t);
% and take LFT with PIE
PIE_est = piess(T,A+Lval*C2,-Lval,{C1(1,:);Kval});
PIE_CL = pielft(PIE,PIE_est);


% =============================================
% === Simulate the system

% % Declare initial values and disturbance
syms st sx real
uinput.ic = [-10*sx; 0];
uinput.w = 10*exp(-st);

% % Set options for discretization and simulation
opts.plot = 'yes';   % don't plot the final solution
opts.N = 8;         % Expand using 8 Chebyshev polynomials
opts.tf = 2;        % Simulate up to t = 2
opts.dt = 1e-2;     % Use time step of 10^-2

% % Perform the actual simulation
% Simulate uncontrolled PIE and extract solution
[solution_OL,grid] = PIESIM(PIE,opts,uinput);
tval = solution_OL.timedep.dtime;
x_OL = reshape(solution_OL.timedep.primary{2}(:,1,:),opts.N+1,[]);
z_OL = solution_OL.timedep.regulated{1}(1,:);
% Simulate controlled PIE and extract solution
[solution_CL,~] = PIESIM(PIE_CL,opts,uinput);
x_CL = reshape(solution_CL.timedep.primary{2}(:,1,:),opts.N+1,[]);
xhat_CL = reshape(solution_CL.timedep.primary{2}(:,2,:),opts.N+1,[]);
z_CL = solution_CL.timedep.regulated{1}(1,:);
zhat_CL = solution_CL.timedep.regulated{1}(3,:);
u_CL = solution_CL.timedep.regulated{1}(2,:);
w = double(subs(uinput.w,st,tval));


echo off


% % Change titles of plots
fig1 = figure(1);  ax1 = fig1.CurrentAxes;
title(ax1,'Open-Loop Observed Output Evolution','Interpreter','latex','FontSize',14);
fig2 = figure(2);  ax2 = fig2.CurrentAxes;
title(ax2,'Open-Loop Regulated Output Evolution','Interpreter','latex','FontSize',14);
subtitle(ax2,'Output $z_1(t)$ and Control Effort $u(t)=z_2(t)$','Interpreter','latex','FontSize',13);
fig3 = figure(3); 
fig3.Children.Title.String = ['Open-Loop Final PDE State $\mathbf{x}_{1}(t=',num2str(opts.tf),',s)$'];
fig4 = figure(4);
fig4.Children(2).Title.String = 'Open-Loop PDE State Evolution $\mathbf{x}_{1}(t,s)$';
fig5 = figure(5);   ax5 = fig5.CurrentAxes;
title(ax5,'Closed-Loop Regulated Output Evolution','Interpreter','latex','FontSize',14);
subtitle(ax5,'True Output $z_1(t)$, Control Effort $u(t)=z_2(t)$, and Estimated Output $z_3(t)$','Interpreter','latex','FontSize',13);
fig6 = figure(6);
fig6.Children(3).String = 'Closed-Loop True ($x_1$) and Estimated ($x_2$) Final PDE State';
fig7 = figure(7);
fig7.Children(5).String = 'Closed-Loop True ($x_1$) and Estimated ($x_2$) PDE State Evolution';

