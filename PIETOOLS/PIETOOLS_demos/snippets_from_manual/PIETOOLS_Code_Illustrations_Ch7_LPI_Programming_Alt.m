% This document illustrates how LPI optimization programs can be declared
% and solved with PIETOOLS. Specifically, we show how an H_infty optimal
% estimator LPI can be solved for a reaction-diffusion equation with
% boundary observations.
% We refer to Chapter 7 of the manual for more context on the codes.
%
% For the example, we consider the following unstable reaction-diffusion
% equation
%  PDEs:     x_{t}(t,s1,s2) = x_{s1s1}(t,s1,s2) +x_{s2s2}(t,s1,s2) + 4x(t,s1,s2) 
%                                   + s1*(1-s1)*(s2+1)*(3-s2)*w(t);
%  With BCs        x(t,0,s2) = 0;      x(t,1,s2) = 0;
%                 x(t,s1,-1) = 0;      x_{s2}(t,s1,1) = 0; 
%  and outputs          z(t) = int_{0}^{1} int_{-1}^{1} x(t,s1,s2) ds1 ds2 + w(t);
%                    y(t,s1) = x(t,s1,1);
% Letting v:=(d^2/ds1^2)(d^2/ds2^2)x, we derive an equivalent PIE of the form:
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
% DB, 12/29/2024: Use pde_var objects instead of sys and state;
% DJ, 01/05/2025: Use 2D example;

clc; clear; close all; clear stateNameGenerator;
%% Declare the system of interest
% Declare the system as a PDE.
pvar s1 s2 t
pde_var state x input w output z sense y;
x.vars = [s1;s2];   x.dom = [0,1;-1,1];
y.vars = s1;        y.dom = [0,1];
PDE = [diff(x,t) == diff(x,s1,2) +diff(x,s2,2) + 4*x + s1*(1-s1)*(s2+1)*(3-s2)*w;
               z == int(x,[s1;s2],[0,1;-1,1]) + w;
               y == subs(x,s2,1);
               subs(x,s1,0)==0;    subs(x,s1,1)==0;
               subs(x,s2,-1)==0;   subs(diff(x,s2),s2,1)==0];

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
prog = lpiprogram([s1;s2],[0,1;-1,1]);
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
Pdeg = 0;
P_opts.exclude = [1;
                  1;1;1;
                  1;1;1;
                  0;1;1;1;1;1;1;1;1];
[prog,P] = poslpivar(prog,Pdim,Pdeg,P_opts);
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
[prog,Z] = lpivar(prog,Zdim,2);


%% 7.3 Imposing Constraints
% % lpi_ineq
% We impose negativity Q<=0 using lpi_ineq:
Iw = eye(size(B1,2));       Iz = eye(size(C1,1));
Q = [-gam*Iw,           -D11',    -(P*B1+Z*D21)'*T;
     -D11,              -gam*Iz,  C1;
     -T'*(P*B1+Z*D21),  C1',      (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];
clear Q_opts
Q_opts.psatz = 1;
prog = lpi_ineq(prog,-Q,Q_opts);

% % lpi_eq
% Alternatively, we could impose negativity Q<=0 manually using lpi_eq:
Rdeg = 3;
[prog_alt,R1] = poslpivar(prog,Q.dim(:,1),Rdeg);
R_opts.psatz = 1;
[prog_alt,R2] = poslpivar(prog_alt,Q.dim(:,1),Rdeg-1,R_opts);
prog_alt = lpi_eq(prog_alt,Q+R1+R2,'symmetric');


%% 7.4 Defining an Objective Function
% Set gam as objective function to MINIMIZE
if isa(gam,'dpvar')
    prog = lpisetobj(prog, gam);
end


%% 7.5 Solving the Optimization Program
% Declare solve options and solve
solve_opts.solver = 'sedumi';         % use SeDuMi to solve the SDP
solve_opts.simplify = true;           % try to simplify the SDP before solving
prog_sol = lpisolve(prog,solve_opts);


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
% Compute the (optimal) controller gain operator K=Z*P^{-1} 
% (not relevant here)
%Kval = getController(Pval,Zval);

% % 7.6.3
% Construct the closed-loop PIE representation of the system using the
% given controller or observer gain.
%PIE_CL = closedLoopPIE(PIE,Kval);
PIE_CL = closedLoopPIE(PIE,Lval,'observer');

