% State-feedback stabilization of an unstable 2D reaction-diffusion PDE.
% Initial coding DB 10/27/2026
% controller reconstruction AT 10/28/2026
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
PIE = convert(PDE_t,'pie');
PIEc = pie2copvar(PIE);
T = PIEc.T; A = PIEc.A; B = PIEc.B2;
vars = {'s1','s2'};
I = eye_copvar_sop(ne,vars,PIE.dom);

%% Positive PI Lyapunov variable and free auxiliary feedback variable.
prog = lpiprogram_sop(T);
[prog,P] = poscopvar(prog,ne,vars,PIE.dom,0);
P = P + 1e-4*I;
spaces.out = []; spaces.in = vars;
[prog,Z] = lpivar_cdopvar(prog,struct('out',ne,'in',ne),spaces,PIE.dom,2);

%% Dual LPI: (A*P+B*Z)*T' + T*(A*P+B*Z)' + 2*alpha*T*P*T' <= 0.
APBZ = A*P+B*Z;
Q = APBZ*T'+T*APBZ'+2*alpha*T*P*T';
prog = lpi_ineq_sop(prog,-Q,struct('like',P,'dD',2,'dw',2));
prog = lpisolve(prog);
info = prog.solinfo.info;
assert(~info.pinf && ~info.dinf && info.numerr==0 && abs(info.feasratio-1)<0.1, ...
    'No control certificate; inspect prog.solinfo.info.');

%% Extract the Lyapunov and auxiliary PI operators.
Psol = getsol_lpivar_sop(prog,P);
Zsol = getsol_lpivar_sop(prog,Z);
fprintf('\nState-feedback LPI feasible for alpha = %.3f; Psol and Zsol extracted.\n',alpha);


%% Controller reconstruction using Quadratic ProgrammingZs
[Ksop,Kop,info] = getController_sop(Psol,Zsol, struct('method', 'qp'));
err = Zsol - Ksop*Psol; norm_err = err*err';

fprintf('\nThe controller is reconstructed. ||Z - KP|| = %.3e\n', full(norm_err.C{1}.params{1}))

 