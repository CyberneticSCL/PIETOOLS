% Periodic reaction-diffusion equation: Jagt et al., Section VI-A.
% Codex, 10/28/2026: Add the constrained stability LPI of Proposition 17.
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
T = PIE.T(2,:); A = PIE.A(2,:); K = PIE.Kop;
%% S*x is the spatial mean; test decay of (I-S)*x.
S = mat2opvar(0,[0,0;1,1],PIE.vars,PIE.dom);
S.R.R1 = 1/2; S.R.R2 = 1/2;
T = (1-S)*T; A = (1-S)*A;
%% The slack X enforces dissipation subject to linear PI constraint
% Initialize lpiprogram
prog = lpiprogram(PIE.vars(:,1),PIE.vars(:,2),PIE.dom);
% Lyapunov decision PI variable (coercive)
[prog,P] = poslpivar(prog,[0;1],3); % Degree at most six.
P = P + 1e-4;
% Lagrange decision PI variable (indefinite)
[prog,X] = lpivar(prog,[1,1;1,0],6);
% Enforce negativity of Lyapunov function D_t<Tx(s,t),PTx(s,t)> <= 0 for
% all x(s,t) : Kx(s,t)=0, a <= s <= b
Q = A'*P*T + T'*P*A + 2*alpha*T'*P*T + X*K + K'*X';
prog = lpi_ineq(prog,-Q,struct('psatz',1));
prog = lpisolve(prog);

info = prog.solinfo.info;
assert(~info.pinf && ~info.dinf && info.numerr==0 && abs(info.feasratio-1)<0.1, ...
    'No stability certificate; inspect prog.solinfo.info.');
Psol = getsol_lpivar(prog,P);
Xsol = getsol_lpivar(prog,X);
fprintf('\nCertified mean-free decay rate alpha = %.4f (analytical limit %.4f).\n', ...
    alpha,pi^2-lam);
