% DEMO10_PBC.m
%
% This document illustrates how the decay rate of a PDE with periodic
% boundary conditions can be certified using PIETOOLS, as in Section VI-A
% of [1].  Specifically, we consider the following system:
% PDE:      d/dt x(t,s) = lam*x(t,s) + (d^2/ds^2) x(t,s),   s in [-1,1];
% BCs:                0 = x(t,-1) - x(t,1);
%                     0 = d/ds x(t,-1) - d/ds x(t,1);
% with lam = pi^2-1. The spatial mean of x grows as exp(lam*t), so the PDE
% is unstable, but the mean-free part x - (1/2)int_{-1}^{1} x(t,s) ds
% decays at rate pi^2-lam = 1.
% First we convert the above PDE to a PIE of the form:
%   [T \dot{xf}](t,s) = [A xf](t,s),    [K xf](t,s) = 0;
% where the PIE state xf includes auxiliary states, and the second row of T
% maps xf to the PDE state x. Letting S be the deviation from the mean,
%   [S x](s) = x(s) - (1/2)int_{-1}^{1} x(r) dr, 
% we then define Tt = S*T(2,:) and At = S*A, so that Tt*xf is the mean-free 
% part of x. To certify that it decays at rate alpha, we solve the LPI
%   find_{P,X}  P,X
%   s.t.        P>=0,
%       At'*P*Tt + Tt'*P*At + 2*alpha*Tt'*P*Tt + X*K + K'*X' <= 0
% where X is a slack variable enforcing the constraint K*xf=0.
%
% The script displays the PDE, solves the LPI for alpha = 0.99 with
% P >= 1e-4*I, stops with an error if the solver reports infeasibility,
% a numerical error or a feasibility ratio more than 0.1 from 1, and prints
% the certified decay rate. Workspace: PDE, PIE, the solved program prog,
% and the solved operators Psol and Xsol.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - DEMO10
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
% DB, 10/29/2026: Initial coding;

% Clear the workspace; declare the variables, the state x on [-1,1], the
% reaction coefficient lam and the decay rate alpha to certify
clc; clear; clear stateNameGenerator
pvar s t;
x = pde_var(s,[-1,1]);
lam = pi^2-1;
alpha = 0.99;                       % Analytical mean-free decay rate: 1.

%% Declare the PDE and convert to the augmented PIE.
eq_PDE = diff(x,t)==lam*x+diff(x,s,2);
eq_BC = [subs(x,s,-1)==subs(x,s,1);
         subs(diff(x,s),s,-1)==subs(diff(x,s),s,1)];
PDE = initialize([eq_PDE;eq_BC]);
display(PDE);
PIE = convert(PDE,'pie');
%% Extract operators defining the PIE
% Declare operator of the mean-free part, 
%   x - (1/2)int_{-1}^{1} x(r) dr, on L2[-1,1]
S = get_op(x - (1/2)*int(x,s,[-1,1]), s);
% Get the PIE dynamics of the mean-free part of the PDE state
T = S*PIE.T(2,:);     A = S*PIE.A(2,:);
% The PIE state must also satisfy K*xf=0
K = PIE.K;
%% The slack X enforces dissipation subject to linear PI constraint
% Initialize lpiprogram
prog = lpiprogram(PIE.vars,PIE.dom);
% Lyapunov decision PI variable (coercive), on the output spaces of T
[prog,P] = poscopvar(prog,T,struct('int',3,'mult',3,'joint',3));
P = P + 1e-4*eye_copvar_sop(T);
% Lagrange decision PI variable (indefinite)
[prog,X] = lpivar_cdopvar(prog,K',6);   % X with the spaces of K'
% Enforce negativity of Lyapunov function D_t<Tx(s,t),PTx(s,t)> <= 0 for
% all x(s,t) : Kx(s,t)=0, a <= s <= b
Q = A'*P*T + T'*P*A + 2*alpha*T'*P*T + X*K + K'*X';
prog = lpi_ineq_sop(prog,-Q);
% Solve the LPI
prog = lpisolve(prog);

% Stop if the solver reports infeasibility or a numerical error
info = prog.solinfo.info;
assert(~info.pinf && ~info.dinf && info.numerr==0 && abs(info.feasratio-1)<0.1, ...
    'No stability certificate; inspect prog.solinfo.info.');
% Extract the solved operators and print the certified decay rate
Psol = getsol_lpivar_sop(prog,P);
Xsol = getsol_lpivar_sop(prog,X);
fprintf('\nCertified mean-free decay rate alpha = %.4f (analytical limit %.4f).\n', ...
    alpha,pi^2-lam);
%[1] D. Jagt, S. Chernyshenko, and M. Peet. Representation and stability analysis of 1D PDEs
% with periodic boundary conditions. In Proceedings of the IEEE Conference on Decision
% and Control, pages 2050–2057. IEEE, 2025.
