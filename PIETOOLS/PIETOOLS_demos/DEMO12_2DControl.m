% DEMO12_2DControl.m
%
% This document illustrates how PIETOOLS can be used for state-feedback
% stabilization of an unstable 2D reaction-diffusion PDE, using 'copvar'
% and 'sopvar' operators.
% Specifically, we consider the following system, with one distributed
% actuator u:
% PDE:      d/dt x(t,s) = lam*x(t,s) + c1*(d^2/ds1^2) x(t,s) + c2*(d^2/ds2^2) x(t,s)
%                           + b(s)*u(t,s),                  s in [0,1]x[0,1];
% BCs:                0 = x(t,s) on the boundary of [0,1]x[0,1] (Dirichlet);
% with c1 = c2 = 1, lam = 2*pi^2+1 and b(s) = s1*(1-s1)*s2*(1-s2).
% The smallest eigenvalue of -(d^2/ds1^2 + d^2/ds2^2) with these BCs is
% 2*pi^2, so the open-loop system is unstable, with growth rate lam-2*pi^2 = 1.
% First we convert the above PDE to a PIE of the form:
%   [T \dot{v}](t) = [A v](t) + [B u](t),
% where v = d^2/(ds1 ds2) x is the fundamental state, and T maps v to x.
% Then, to find a state feedback u = K*T*v that makes the closed-loop system
% decay at rate alpha, we solve the dual LPI
%   find_{P,Z}  P,Z
%   s.t.        P>=0,
%       (A*P+B*Z)*T' + T*(A*P+B*Z)' + 2*alpha*T*P*T' <= 0
% where P is a positive operator, and Z is a free operator declared using
% `lpivar_cdopvar`. A feedback law K is then reconstructed from Z = K*P
% using `getController_sop`.
%
% The PIE that `convert` returns is a 'pie_struct' whose operators are
% 'copvar' objects. The LPI program is declared on PIE.vars and
% PIE.dom; P takes its spaces and domain from T', on which it acts in
% T*P*T', and Z those of B', since Z = K*P maps the state to the input.
%
% The script solves the LPI for alpha = 0.1 with P >= 1e-4*I, stops with
% an error if the solver reports infeasibility, a numerical error or a
% feasibility ratio more than 0.1 from 1, reconstructs K by quadratic
% programming, and prints a residual of Z = K*P. Workspace: PIE, the
% solved program prog, Psol, Zsol, and the controller Ksop, Kop.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO12
%
% Copyright (C)2026  PIETOOLS Team
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
% DB 10/27/2026 Initial coding 

% Clear the workspace; declare the variables, the coefficients, the
% number of state components ne and the decay rate alpha
clear; clc; clear stateNameGenerator
pvar s1 s2
c1 = 1; c2 = 1; lam = 2*pi^2+1; ne = 1;
alpha = 0.1;

%% PDE with one distributed actuator per state component.
x = pde_var(ne,[s1;s2],[0,1;0,1]);
u = pde_var('control',ne);
b = s1*(1-s1)*s2*(1-s2);       
PDE_t = [diff(x,'t')==lam*x+c1*diff(x,s1,2)+c2*diff(x,s2,2)+b*u;
         subs(x,s1,0)==0; subs(x,s1,1)==0;
         subs(x,s2,0)==0; subs(x,s2,1)==0];
% Convert to a PIE with 'copvar' operators; extract T, A and Bu
PIE = convert(PDE_t,'pie');
T = PIE.T; A = PIE.A; B = PIE.Bu;

%% Positive PI Lyapunov variable and free auxiliary feedback variable.
prog = lpiprogram(PIE.vars,PIE.dom);
[prog,P] = poscopvar(prog,T',0);        % P on the input spaces of T
P = P + 1e-4*eye_copvar_sop(P);         % enforce P>=1e-4
[prog,Z] = lpivar_cdopvar(prog,B',2);   % Z with the spaces of B'

%% Dual LPI: (A*P+B*Z)*T' + T*(A*P+B*Z)' + 2*alpha*T*P*T' <= 0.
APBZ = A*P+B*Z;
Q = APBZ*T'+T*APBZ'+2*alpha*T*P*T';
prog = lpi_ineq_sop(prog,-Q,struct('dD',2,'dw',2));  % like: P, inferred
% Solve; stop if the solver reports infeasibility or a numerical error
prog = lpisolve(prog);
info = prog.solinfo.info;
assert(~info.pinf && ~info.dinf && info.numerr==0 && abs(info.feasratio-1)<0.1, ...
    'No control certificate; inspect prog.solinfo.info.');

%% Extract the Lyapunov and auxiliary PI operators.
Psol = getsol_lpivar_sop(prog,P);
Zsol = getsol_lpivar_sop(prog,Z);
fprintf('\nState-feedback LPI feasible for alpha = %.3f; Psol and Zsol extracted.\n',alpha);


%% Controller reconstruction using quadratic programming
[Ksop,Kop,info] = getController_sop(Psol,Zsol, struct('method', 'qp'));
% Residual err = Z - K*P, a map L2 -> R^ne: ||err||^2 is the 2-norm of the
% ne x ne matrix err*err'
err = Zsol - Ksop*Psol; norm_err = err*err';

fprintf('\nThe controller is reconstructed. ||Z - KP|| = %.3e\n', sqrt(norm(full(norm_err.C{1}.params{1}))))

 