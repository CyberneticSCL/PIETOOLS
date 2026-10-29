% Periodic reaction-diffusion equation from Jagt et al., 2025 Section VI-A.
clc; clear; clear stateNameGenerator
pvar s t;
x = pde_var(s,[-1,1]);
lam = pi^2-1;
alpha = 0.99;                       % Analytical mean-free decay rate: 1.
% Initial coding DB 10/29/2026

%% Declare the PDE and convert to the augmented PIE.
eq_PDE = diff(x,t)==lam*x+diff(x,s,2);
eq_BC = [subs(x,s,-1)==subs(x,s,1);
         subs(diff(x,s),s,-1)==subs(diff(x,s),s,1)];
PDE = initialize([eq_PDE;eq_BC]);
display(PDE);
PIE = convert(PDE,'pie');
display(PIE);
PIEc = pie2copvar(PIE);
%% convert to copvar
% First row corresponds to te average/moment state, second row corresponds to the original PDE state
T = copvar(PIEc.T.C(2,:)); A = copvar(PIEc.A.C(2,:));
F = PIEc.F; K = PIEc.K;
%% S*x is the spatial mean; test decay of (I-S)*x.
S = copvar(PIEc.T.C(2,1))*F;
I = eye_copvar_sop(1,{'s'},PIE.dom);
T = (I-S)*T; A = (I-S)*A;
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
% D. S. Jagt and S. Chernyshenko and M. M. Peet, "Representation and Stability Analysis of 1D PDEs with Periodic
% Boundary Conditions," arXiv:2503.22896, 2025.