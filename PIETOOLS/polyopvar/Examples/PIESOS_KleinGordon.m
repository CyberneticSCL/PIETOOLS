%
% Kelin-Gordon equation with stabilizing nonlinearity
%
% PDE: u_{tt} = u_{xx} -u_{t}*u_{x}^2;
% BCs: u(t,0) = u_{x}(t,1) = 0;
%
% Introduce u1 = u_{t} and u2 = u_{x}, then
%       u1_{t} = u2_{x} -u1*u2^2;
%       u2_{t} = u1_{x};
%       u1(t,0) = 0;    u2(t,1) = 0;
%

clear;  clear stateNameGenerator
pvar x t

% Declare the nonlinear PDE
u1 = pde_var(x,[0,1]);
u2 = pde_var(x,[0,1]);
PDE = [diff(u1,t)==diff(u2,x)-u1*u2^2;
       diff(u2,t)==diff(u1,x);
       subs(u1,x,0)==0;  subs(u2,x,1)==0];

% Convert to a PIE
PIE = convert(PDE);
T = PIE.T;  f = PIE.f;
u = f.vartab;

% Initialize a distributed SOS program
prog = piesos_program(u);

% Declare a distributed SOS candidate Lyapunov functional
Zd = dmonomials(u,1);
[prog,V] = piesos_sosvar(prog,Zd);
V = V + 1e-2*innerprod(u,u);
                                              
% Enforce dV <= -s*g <= 0 along solutions of the PIE
dV = Liediff(V,PIE); 

%ineq_opts.psatz = 0:1;
prog = piesos_ineq(prog,-dV);

% Solve the distributed sos program
prog = piesos_solve(prog);