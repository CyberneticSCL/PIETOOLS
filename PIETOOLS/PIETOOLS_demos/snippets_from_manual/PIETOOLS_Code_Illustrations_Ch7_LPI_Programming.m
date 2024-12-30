
% This document illustrates how LPI optimization programs can be declared
% and solved with PIETOOLS. Specifically, we show how an H_infty optimal
% estimator LPI can be solved for a reaction-diffusion equation with
% boundary observations.
% We refer to Chapter 7 of the manual for more context on the codes.
%
% For the example, we consider the following unstable reaction-diffusion
% equation
%  PDEs:            x_{t}(t) = x_{ss}(t,s) + 4x(t,s) + w(t);
%  With BCs           x(t,0) = 0;      x_{s}(t,1) = 0;
%  and outputs          z(t) = int_{0}^{1} x(t,s) ds + w(t);
%                       y(t) = x(t,1);
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
% We show how this LPI can be solved here.
% DB, 12/29/2024: Use pde_var objects instead of sys and state

clc; clear; close all; clear stateNameGenerator;
%% Declare the system of interest
% Declare the system as a PDE.
pvar s t
x=pde_var('state',1,s,[0 1]);
y=pde_var('sense',1); z=pde_var('output',1);
w=pde_var('input',1);
PDE = [diff(x,t) == diff(x,s,2) + 4*x + w;
             z == int(x,s,[0,1]) + w;
             y == subs(x,s,1);
             subs(x,s,0) == 0;
             subs(diff(x,s),s,1) == 0];
% Convert to PIE and extract the relevant operators.
PIE = convert(PDE,'pie');
T = PIE.T;
A = PIE.A;    C1 = PIE.C1;    C2 = PIE.C2;
B1 = PIE.B1;  D11 = PIE.D11;  D21 = PIE.D21;


%% 7.1 Initializing an Optimization Program
% An LPI optimization program must be initialized using "lpiprogram". This
% function requires the independent variables and their domain to be
% specified, and also allows decision variables to be declared (though
% those can be added later as well)
%pvar s1 s2;                            % Declare free (polynomial) variables
%dom = [0,1;-1,1];                      % Declare domain 
%dpvar d1 d2;                           % Initialize decision variables
%prog = lpiprogram([s1;s2], dom, [d1;d2]);  
% In our case
prog = lpiprogram(s,[0,1]);
prog.vartable
PIE.vars


%% 7.2 Declaring Decision Variables
% % 7.2.1 lpidecvar
% Decision variables with a desired name can be declared using "dpvar"
%dpvar d1;                              % Initialize a new decision variable
%prog1 = lpidecvar(prog,d1);            % Add the decision variable to the program
% or directly using "lpidecvar"
%[prog2,d1] = lpidecvar(prog,'d1');     % Generate decision variable with name d1
% "lpidecvar" can also be used to declare an array of decision variables
%m = 2;  n = 3;
%[prog3,d_arr] = lpidecvar(prog,[m,n]); % Generate mxn array of decision variables
% In our program, we only have 1 scalar decision variable:
[prog,gam] = lpidecvar(prog,'gam');
% Alternatively, we can also construct an observer for given gam
%gam = 2;

% % 7.2.2 poslpivar
% "poslpivar(prog,Pdim,Pdeg,opts)" generates a positive semidefintie PI
% operator decision variable of dimension Pdim x Pdim, using monomials of
% degrees specified by Pdeg, and using additional options opts. Only the
% field "Pdim" is mandatory.
Pdim = T.dim(:,1);
Pdeg = {4,[2,3,4],[2,3,4]};
opts.sep = 1;
[prog,P] = poslpivar(prog,Pdim,Pdeg,opts);
% NOTE: "poslpivar" only generates positive semidefinite operators. To
% ensure strict positivity, add a positive operator eppos*I:
eppos = 1e-4;
P = P +eppos*eye(size(P));

% % 7.2.3 lpivar
% "lpivar(prog,Zdim,Zdeg)" generates a PI operator decision variable 
% of dimension Zdim, using monomials of
% degrees specified by Pdeg, and using additional options opts. Only the
% field "Pdim" is mandatory.
Zdim = C2.dim(:,[2,1]);
Zdeg = [4,0,0];
[prog,Z] = lpivar(prog,Zdim,Zdeg);


%% 7.3 Imposing Constraints
% % lpi_ineq
% We impose negativity Q<=0 using lpi_ineq:
Iw = eye(size(B1,2));       Iz = eye(size(C1,1));
Q = [-gam*Iw,           -D11',    -(P*B1+Z*D21)'*T;
     -D11,              -gam*Iz,  C1;
     -T'*(P*B1+Z*D21),  C1',      (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];
prog = lpi_ineq(prog,-Q);

% % lpi_eq
% Alternatively, we could impose negativity Q<=0 manually using lpi_eq:
Rdeg =  {6,[4,5,6],[4,5,6]};
[prog_alt,R] = poslpivar(prog,Q.dim(:,1),Rdeg);
prog_alt = lpi_eq(prog_alt,Q+R,'symmetric');


%% 7.4 Defining an Objective Function
% Set gam as objective function to MINIMIZE
if isa(gam,'dpvar')
    prog = lpisetobj(prog, gam);
end


%% 7.5 Solving the Optimization Program
% Declare solve options and solve
opts.solver = 'sedumi';         % use SeDuMi to solve the SDP
opts.simplify = true;           % try to simplify the SDP before solving
prog_sol = lpisolve(prog,opts);


%% 7.6 Extracting the Solution
% Extract solved (optimal) values of the decision variable gam and decision
% variable operators P and Z;
gam_val = lpigetsol(prog_sol,gam)
Pval = lpigetsol(prog_sol,P);
Zval = lpigetsol(prog_sol,Z);

% % 7.6.1 getObserver
% Compute the (optimal) observer gain operator L=P^{-1}*Z
Lval = getObserver(Pval,Zval);

% % 7.6.2 getController
% Compute the (optimal) controler gain operator K=Z*P^{-1} 
% (not relevant here)
%Kval = getController(Pval,Zval);

% % 7.6.3
% Construct the closed-loop PIE representation of the system using the
% given controller or observer gain.
%PIE_CL = closedLoopPIE(PIE,Kval);
PIE_CL = closedLoopPIE(PIE,Lval,'observer');

