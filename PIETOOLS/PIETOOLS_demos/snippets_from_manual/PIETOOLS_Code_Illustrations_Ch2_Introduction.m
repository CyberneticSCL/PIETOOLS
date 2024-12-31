% This document illustrates how PIETOOLS can be used to simulate system
% dynamics, analyse stability and design optimal controllers.
% We refer to Chapter 2 of the manual for more context on the codes.
%  
% For this example, we consider the damped wave equation with dynamic 
% boundary input,
%                xi_{tt}=c xi_{ss}-b xi_{t}+sw(t),
% and with regulated output
%               r(t) = int_0^1 xi_{s} ds = xi(s=1)-xi(s=0).
% Introducing the state phi=(phi1,phi2) with phi1=xi_{s} and phi2=xi_{t},
% this wave equation can be expressed as an ODE-PDE system
%  ODE                 x_{t} = -x + u(t) 
%  PDEs:            phi1_{t} = phi2_{s}
%                   phi2_{t} = c phi1_{s} - b phi2_{t}+sw(t)
%  With BCs        phi1(s=0) = 0
%                  phi2(s=1) = x = exp(-t)x(t=0)+ int_{0}^{t} exp(-t+tau) u(tau) d tau
%  and outputs         z1(t) = int_{0}^{1} phi1(t,s) ds
%                      z2(t) = u(t)
% where now we add an output z2(t)=u(t) so that, in constructing an optimal
% control, we also minimize the control effort.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Code Illustrations
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
% If you modify this code, make sure to change the code in the manual as
% well, and vice versa. Document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 12/28/2024: Initial coding;
% DB, 12/29/2024: Use pde_var objects instead of sys and state
clear


%% 2.2.1 Declare the ODE-PDE model
% Declare an ODE state variable x, a vector-valued 1D PDE state variable
% phi with 2 elements, and a finite-dimensional exogenous input w,
% controlled input u, and regulated output z, with z containing 2 elements. 
clc; clear; close all; clear stateNameGenerator;
% Declare independent variables
pvar t s;   
% Declare state, input, and output variables
phi = pde_var(2,s,[0,1]);   x = pde_var();
w = pde_var('in');          u = pde_var('control');   
z = pde_var('out',2);
% Declare system parameters
c=1;    b=.01;
% Declare ODE-PDE equations
%   d/dt      x(t) = -x(t) + u(t);
%   d/dt phi1(t,s) = phi2_{ss}(t,s);
%   d/dt phi2(t,s) = c*phi1_{ss}(t,s) + s*w(t) -b*phi2(t,s);
eq_dyn = [diff(x,t,1)==-x+u
          diff(phi,t,1)==[0 1; c 0]*diff(phi,s,1)+[0;s]*w+[0 0;0 -b]*phi];
% Declare output equations
%            z1(t) = int_{0}^{1} phi1(t,s)ds;
%            z2(t) = u(t);
eq_out= z==[int([1 0]*phi,s,[0,1]); 
            u];

%% 2.2.2 Declare the boundary conditions
bc1 = [0 1]*subs(phi,s,0) == 0;         % phi2(t,0) = 0;
bc2 = [1 0]*subs(phi,s,1) == x;         % phi1(t,1) = x(t);
% create the PDE system
odepde = [eq_dyn;eq_out;bc1;bc2];
odepde = initialize(odepde);

%% 2.2.3: Simulate the ODE-PDE system

% Declare simulation settings.
opts.plot = 'no';   % don't plot final solution
opts.N = 8;         % expand solution using 8 Chebyshev polynomials
opts.tf = 12;       % simulate up to t = 15
opts.dt = 3*1e-2;   % use time step of 3*10^-2
opts.intScheme = 1; % use backward differentiation formula for time integration

% Declare initial values and values of disturbance and input.
syms st sx real                 % symbolic time and space variables
uinput.ic.PDE = [0,0];          % initial values of both PDE state variables
uinput.ic.ODE = 0;              % initial value of ODE state variable
uinput.u = 0;                   % value of control as function of time
uinput.w = sin(5*st)*exp(-st);  % value of disturbance as function of time

% Declare number of PDE state variables differentiable up to each order:
%   0 0th-order, 2 1st-order, and 0 2nd-order differentiable states.
ndiff = [0,2,0]; 

% Simulate the ODE-PDE system using specified settings.
[solution,grids] = PIESIM(odepde, opts, uinput, ndiff);

% Extract values of the second state variables and thet output at each time step.
tval = solution.timedep.dtime;
phi2 = reshape(solution.timedep.pde(:,2,:),opts.N+1,[]);
zval = solution.timedep.regulated;
wval = subs(uinput.w,st,tval);

% Plot evolution of the PDE state and regulated output.
fig1 = figure('Position',[200 150 1000 450]);
set(gcf, 'Color', 'w');
sgtitle('Open loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',16);

ax1 = subplot(1,2,1);
box on
surf(tval,grids.phys,phi2,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar;     
colormap jet
ylabel(h,'$|\dot{\mathbf{x}}(t,s)|$','interpreter', 'latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
title('PDE state $\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
ax1.Position = [0.08,0.12,0.33,0.73];
ax1.TickLabelInterpreter = 'latex';
h.Label.String = '';        h.TickLabelInterpreter = 'latex';
ax1.ZLim = [-0.21,0.21];

ax2 = subplot(1,2,2);
plot(tval,wval,'k',tval,zval(1,:),'r','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$w(t)$','$r(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$r(t)$','FontSize',15,'Interpreter','latex');
title('Regulated output $r(t)$','Interpreter','latex','FontSize',15);
ax2.Position = [0.63,0.12,0.33,0.73];
ax2.TickLabelInterpreter = 'latex';
ax2.YLim = [-0.6,0.8];
%saveas(fig1,'Ch2_OpenLoop_Plot','epsc');


%% 2.2.4: Analyse stability and construct a controller

% Convert the ODE-PDE to a PIE.
PIE = convert(odepde,'pie');

% Verify stability of the system by runnig the stability executive.
settings = lpisettings('heavy');    % use heavy settings: more expensive but less conservative
lpiscript(PIE,'stability',settings);

% Compute an upper bound gam on the H_infty norm or L2-gain of the system,
% so that ||z(.)||_{L2} <= gam*||w(.)||_{L2}.
[~,~, gam] = lpiscript(PIE,'l2gain',settings);

% Compute a feedback gain operator K that minimizes the H_infty norm of the
% closed-loop system with state feedback u = K*v, for fundamental (PIE)
% state v.
[~, Kval, gam_CL] = lpiscript(PIE,'hinf-controller','light');

% Now, we construct the closed-loop system using the controller obtained by
% solving the above LPI and then re-runs the simulations to see if there
% was any improvement in the performance.
PIE_CL = closedLoopPIE(PIE,Kval);
PIE_CL = pie_struct(PIE_CL);
PIE_CL = initialize(PIE_CL);


%% -.-.- Simulate the closed-loop ODE-PDE system
% Simulate the solution to the PIE with controller
[solution_CL,grids_CL] = PIESIM(PIE_CL,opts,uinput,ndiff);

% Having found an optimal controller, constructed closed loop system, and
% running PIESIM to simulate, we can now extract the solution and plot it
% against the open-loop response to see the benefit in using the
% controller. Firstly, we notice that the L2-gain bound has significantly
% lowered (from 5.2 to 0.78). Then, looking at the output in presence of a
% disturbance, we see that the effect of disturbance on this neutrally
% stable system is reduced.
tval = solution_CL.timedep.dtime;
phi1_cl = reshape(solution_CL.timedep.pde(:,1,:),opts.N+1,[]);
phi2_cl = reshape(solution_CL.timedep.pde(:,2,:),opts.N+1,[]);
zval_cl = solution_CL.timedep.regulated;
wval = subs(uinput.w,st,tval);



% Plot closed-loop evolution of the PDE state and regulated output.
fig2 = figure('Position',[200 150 1000 450]);
set(gcf, 'Color', 'w');
sgtitle('Closed loop zero-state response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',16);

ax1 = subplot(1,2,1);
box on
surf(tval,grids.phys,phi2_cl,'FaceAlpha',0.75,'Linestyle','--','FaceColor','interp','MeshStyle','row');
h=colorbar;     
colormap jet
ylabel(h,'$|\dot{\mathbf{x}}(t,s)|$','interpreter', 'latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    ylabel('$s$','FontSize',15,'Interpreter','latex');
zlabel('$\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
title('PDE state $\dot{\mathbf{x}}(t,s)$','FontSize',15,'Interpreter','latex');
ax1.Position = [0.08,0.12,0.33,0.73];
ax1.TickLabelInterpreter = 'latex';
h.Label.String = '';    h.TickLabelInterpreter = 'latex';
ax1.ZLim = [-0.21,0.21];

ax2 = subplot(1,2,2);
plot(tval,wval,'k',tval,zval_cl(1,:),'r',tval,zval_cl(2,:),'b','LineWidth',2)
grid on
box on
set(gcf, 'Color', 'w');
legend('$w(t)$','$r(t)$','$u(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$r(t)$','FontSize',15,'Interpreter','latex');
title('Regulated output $r(t)$','Interpreter','latex','FontSize',15);
ax2.Position = [0.63,0.12,0.33,0.73];
ax2.TickLabelInterpreter = 'latex';
ax2.YLim = [-0.6,0.8];
%saveas(fig2,'Ch2_ClosedLoop_Plot','epsc');


% Compare evolution of open-loop and closed-loop regulated output r(t).
fig3 = figure('Position',[200 150 1000 350]);
set(gcf, 'Color', 'w');
box on
plot(tval,zval(1,:),'r--',tval,zval_cl(1,:),'r','LineWidth',2)
legend('$r(t)$','$r_{cl}(t)$','Interpreter','latex','FontSize',15)
xlabel('$t$','FontSize',15,'Interpreter','latex');    
ylabel('$r(t)$','FontSize',15,'Interpreter','latex');
title('Open- and closed-loop zero-state regulated output response with $w(t)=sin(5t)e^{-t}$','Interpreter','latex','FontSize',15);
grid on
set(gca,'TickLabelInterpreter','latex');
%saveas(fig3,'Ch2_RegulatedOutput_Plot','epsc');