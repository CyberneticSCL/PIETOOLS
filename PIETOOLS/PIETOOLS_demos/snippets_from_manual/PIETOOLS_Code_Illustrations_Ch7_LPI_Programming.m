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
% T = PIE.T;                                                               % DB, 11/04/2026 (was)
% A = PIE.A;    C1 = PIE.C1;    C2 = PIE.C2;                                % DB, 11/04/2026 (was)
% B1 = PIE.B1;  D11 = PIE.D11;  D21 = PIE.D21;                              % DB, 11/04/2026 (was)
% PIEo = pie2opvar(PIE);                                                    % DB, 11/04/2026
% T = PIEo.T;                                                               % DB, 11/04/2026
% A = PIEo.A;    C1 = PIEo.C1;    C2 = PIEo.C2;                             % DB, 11/04/2026
% B1 = PIEo.B1;  D11 = PIEo.D11;  D21 = PIEo.D21;                           % DB, 11/04/2026
% The operators of the PIE are 'copvar' objects, used directly below       % DB, 11/06/2026
T = PIE.T;                                                                  % DB, 11/06/2026
A = PIE.A;    C1 = PIE.C1;    C2 = PIE.C2;                                  % DB, 11/06/2026
B1 = PIE.B1;  D11 = PIE.D11;  D21 = PIE.D21;                                % DB, 11/06/2026
dom = struct('vars',{{'s1','s2'}},'dom',[0,1;-1,1]);                        % DB, 11/06/2026


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
% Pdim = T.dim(:,1);                                                        % DB, 11/06/2026 (was)
% Pdeg = 0;                                                                 % DB, 11/06/2026 (was)
% P_opts.exclude = [1;                                                      % DB, 11/06/2026 (was)
%                   1;1;1;                                                  % DB, 11/06/2026 (was)
%                   1;1;1;                                                  % DB, 11/06/2026 (was)
%                   0;1;1;1;1;1;1;1;1];                                     % DB, 11/06/2026 (was)
% [prog,P] = poslpivar(prog,Pdim,Pdeg,P_opts);                              % DB, 11/06/2026 (was)
% NOTE: "poslpivar" only generates positive semidefinite operators. To
% ensure strict positivity, add a positive operator eppos*I:
% eppos = 1e-4;                                                             % DB, 11/06/2026 (was)
% P = P +eppos*eye(size(P));                                                % DB, 11/06/2026 (was)
% For 'copvar' objects, the spaces and dimensions of the state are read     % DB, 11/06/2026
% from T, and "poscopvar(prog,Pdims,Pspaces,dom,Pdeg,opts)" declares P      % DB, 11/06/2026
% on them. Here, we keep only the multiplier term (multi-index [1 1]).      % DB, 11/06/2026
[Pspaces,Pdims] = copvar_space_list(T,'out');                               % DB, 11/06/2026
Pdeg = 0;                                                                   % DB, 11/06/2026
P_opts.include = {[1 1]};                                                   % DB, 11/06/2026
[prog,P] = poscopvar(prog,Pdims,Pspaces,dom,Pdeg,P_opts);                   % DB, 11/06/2026
eppos = 1e-4;                                                               % DB, 11/06/2026
P = P + eppos*eye_copvar_sop(Pdims,Pspaces,dom);                            % DB, 11/06/2026

% % 7.2.3 lpivar
% "lpivar(prog,Zdim,Zdeg)" generates a PI operator decision variable 
% of dimension Zdim, using monomials of
% degrees specified by Pdeg, and using additional options opts. Only the
% field "Pdim" is mandatory.
% Zdim = C2.dim(:,[2,1]);                                                   % DB, 11/06/2026 (was)
% Zdeg = 2;                                                                 % DB, 11/06/2026 (was)
% [prog,Z] = lpivar(prog,Zdim,Zdeg);                                        % DB, 11/06/2026 (was)
% For 'copvar' objects, "lpivar_cdopvar" takes the output and input spaces  % DB, 11/06/2026
% and dimensions of Z (here those of the input and output of C2).           % DB, 11/06/2026
[Zout_sp,Zout_dm] = copvar_space_list(C2,'in');                             % DB, 11/06/2026
[Zin_sp,Zin_dm] = copvar_space_list(C2,'out');                              % DB, 11/06/2026
Zdeg = 2;                                                                   % DB, 11/06/2026
[prog,Z] = lpivar_cdopvar(prog,struct('out',Zout_dm,'in',Zin_dm), ...       % DB, 11/06/2026
                          struct('out',{Zout_sp},'in',{Zin_sp}),dom,Zdeg);  % DB, 11/06/2026


%% 7.3 Imposing Constraints
% % lpi_ineq
% We impose negativity Q<=0 using lpi_ineq:
% Iw = eye(size(B1,2));       Iz = eye(size(C1,1));                         % DB, 11/06/2026 (was)
% Q = [-gam*Iw,           -D11',    -(P*B1+Z*D21)'*T;                       % DB, 11/06/2026 (was)
%      -D11,              -gam*Iz,  C1;                                     % DB, 11/06/2026 (was)
%      -T'*(P*B1+Z*D21),  C1',      (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];           % DB, 11/06/2026 (was)
% clear Q_opts                                                              % DB, 11/06/2026 (was)
% Q_opts.psatz = 1;                                                         % DB, 11/06/2026 (was)
% prog = lpi_ineq(prog,-Q,Q_opts);                                          % DB, 11/06/2026 (was)
[Iw_sp,Iw_dm] = copvar_space_list(B1,'in');                                 % DB, 11/06/2026
[Iz_sp,Iz_dm] = copvar_space_list(C1,'out');                                % DB, 11/06/2026
Iw = eye_copvar_sop(Iw_dm,Iw_sp,dom);                                       % DB, 11/06/2026
Iz = eye_copvar_sop(Iz_dm,Iz_sp,dom);                                       % DB, 11/06/2026
D21r = on_registry_sop(D21,Z);      % D21 on the registry of Z             % DB, 11/06/2026
Q = [-gam*Iw,           -D11',    -(P*B1+Z*D21r)'*T;                        % DB, 11/06/2026
     -D11,              -gam*Iz,  C1;                                       % DB, 11/06/2026
     -T'*(P*B1+Z*D21r), C1',      (P*A+Z*C2)'*T+T'*(P*A+Z*C2)];             % DB, 11/06/2026
% "lpi_ineq_sop" declares a positive slack, with monomial degrees read from % DB, 11/06/2026
% the support of Q (lift degree and weight degree raised by dD, dw = 1)      % DB, 11/06/2026
clear Q_opts                                                                % DB, 11/06/2026
Q_opts.like = P;    Q_opts.dD = 2;                                          % DB, 11/06/2026
prog = lpi_ineq_sop(prog,-Q,Q_opts);                                        % DB, 11/06/2026

% % lpi_eq
% Alternatively, we could impose negativity Q<=0 manually using lpi_eq:
% Rdeg = 3;                                                                 % DB, 11/06/2026 (was)
% [prog_alt,R1] = poslpivar(prog,Q.dim(:,1),Rdeg);                          % DB, 11/06/2026 (was)
% R_opts.psatz = 1;                                                         % DB, 11/06/2026 (was)
% [prog_alt,R2] = poslpivar(prog_alt,Q.dim(:,1),Rdeg-1,R_opts);             % DB, 11/06/2026 (was)
% prog_alt = lpi_eq(prog_alt,Q+R1+R2,'symmetric');                          % DB, 11/06/2026 (was)
Rdeg = 3;                                                                   % DB, 11/06/2026
[Q_sp,Q_dm] = copvar_space_list(Q,'out');                                   % DB, 11/06/2026
[prog_alt,R1] = poscopvar(prog,Q_dm,Q_sp,dom,Rdeg);                         % DB, 11/06/2026
R_opts.psatz = 1;                                                           % DB, 11/06/2026
[prog_alt,R2] = poscopvar(prog_alt,Q_dm,Q_sp,dom,Rdeg-1,R_opts);            % DB, 11/06/2026
prog_alt = lpi_eq_sop(prog_alt,Q+R1+R2,'symmetric');                        % DB, 11/06/2026


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
% gam_val = lpigetsol(prog_sol,gam)                                         % DB, 11/06/2026 (was)
% Pval = lpigetsol(prog_sol,P);                                             % DB, 11/06/2026 (was)
% Zval = lpigetsol(prog_sol,Z);                                             % DB, 11/06/2026 (was)
gam_val = double(lpigetsol_sop(prog_sol,gam))                               % DB, 11/06/2026
Pval = getsol_lpivar_sop(prog_sol,P);                                       % DB, 11/06/2026
Zval = getsol_lpivar_sop(prog_sol,Z);                                       % DB, 11/06/2026

% % 7.6.1 getObserver
% Compute the (optimal) observer gain operator L=P^{-1}*Z
% Lval = getObserver(Pval,Zval);                                            % DB, 11/06/2026 (was)
% P is separable (a multiplier plus eppos*I), so P^{-1}*Z is computed exactly % DB, 11/06/2026
Lval = getObserver_sop(Pval,Zval);                                          % DB, 11/06/2026

% % 7.6.2 getController
% Compute the (optimal) controller gain operator K=Z*P^{-1} 
% (not relevant here)
%Kval = getController(Pval,Zval);

% % 7.6.3
% Construct the closed-loop PIE representation of the system using the
% given controller or observer gain.
%PIE_CL = closedLoopPIE(PIE,Kval);
PIE_CL = closedLoopPIE(PIE,Lval,'observer');

