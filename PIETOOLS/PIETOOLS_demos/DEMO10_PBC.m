% DEMO10_PBC.m
%
% This document illustrates how the decay rate of a PDE with periodic
% boundary conditions can be certified using PIETOOLS, as in Section VI-A
% of D. S. Jagt, S. Chernyshenko and M. M. Peet, "Representation and
% Stability Analysis of 1D PDEs with Periodic Boundary Conditions,"
% arXiv:2503.22896, 2025.
% Specifically, we consider the following system:
% PDE:      d/dt x(t,s) = lam*x(t,s) + (d^2/ds^2) x(t,s),   s in [-1,1];
% BCs:                0 = x(t,-1) - x(t,1);
%                     0 = d/ds x(t,-1) - d/ds x(t,1);
% with lam = pi^2-1. The spatial mean of x grows as exp(lam*t), so the PDE
% is unstable, but the mean-free part x - (1/2)int_{-1}^{1} x(t,s) ds
% decays at rate pi^2-lam = 1.
% First we convert the above PDE to a PIE of the form:
%   [T \dot{v}](t,s) = [A v](t,s),    [K v](t,s) = 0;
% where the PIE state v includes auxiliary states, and the second row of T
% maps v to the PDE state x. Letting S be the mean operator,
% [S x](s) = (1/2)int_{-1}^{1} x(r) dr, we then define Tt = (I-S)*T and
% At = (I-S)*A, so that Tt*v is the mean-free part of x. To certify that
% it decays at rate alpha, we solve the LPI
%   find_{P,X}  P,X
%   s.t.        P>=0,
%       At'*P*Tt + Tt'*P*At + 2*alpha*Tt'*P*Tt + X*K + K'*X' <= 0
% where X is a slack variable enforcing the constraint K*v=0.
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
% DB, 10/29/2026: Initial coding;
% DJ, 10/30/2026: Update to use "int2copvar" to compute mean;

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
display(PIE);
%% Extract operators defining the PIE
% First row of T corresponds to auxiliary state variables,
% second row corresponds to the original PDE state
T = copvar(PIE.T.C(2,:));   A = copvar(PIE.A.C(2,:));
% The PIE state must also satisfy K*x=0
K = PIE.K;
%% Take the deviation of the PDE state from the mean
S = (1/2)*int2copvar(eye(1),T.vars,T.dom);  % S*x = 0.5*int_{-1}^{1} x(s) ds
I = eye_copvar_sop(1,T.vars,T.dom);        % I*x = x;
T = (I-S)*T;            A = (I-S)*A;
% %% S*x is the spatial mean; test decay of (I-S)*x.
% S = copvar(PIE.T.C(2,1))*F;
% I = eye_copvar_sop(1,{'s'},PIE.dom);
% T = (I-S)*T; A = (I-S)*A;
%% The slack X enforces dissipation subject to linear PI constraint
% Initialize lpiprogram
prog = lpiprogram_sop(T);
% Lyapunov decision PI variable (coercive)
[prog,P] = poscopvar(prog,1,{'s'},PIE.dom,struct('int',3,'mult',3,'joint',3));
P = P + 1e-4*I;
% Lagrange decision PI variable (indefinite)
spaces.out = {[], {'s'}};  % R and L2[s]
spaces.in  = [];         % R
[prog,X] = lpivar_cdopvar(prog,struct('out',[1;1],'in',1),spaces,PIE.dom,6);
% Enforce negativity of Lyapunov function D_t<Tx(s,t),PTx(s,t)> <= 0 for
% all x(s,t) : Kx(s,t)=0, a <= s <= b
Q = A'*P*T + T'*P*A + 2*alpha*T'*P*T + X*K + K'*X';
prog = lpi_ineq_sop(prog,-Q);
prog = lpisolve(prog);

info = prog.solinfo.info;
assert(~info.pinf && ~info.dinf && info.numerr==0 && abs(info.feasratio-1)<0.1, ...
    'No stability certificate; inspect prog.solinfo.info.');
Psol = getsol_lpivar_sop(prog,P);
Xsol = getsol_lpivar_sop(prog,X);
fprintf('\nCertified mean-free decay rate alpha = %.4f (analytical limit %.4f).\n', ...
    alpha,pi^2-lam);