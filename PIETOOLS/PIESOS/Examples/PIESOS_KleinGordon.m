%
% Klein-Gordon equation with stabilizing nonlinearity
%
% PDE: u_{tt} = u_{xx} - u_{t}*u_{x}^2;
% BCs: u(t,0) = u_{x}(t,1) = 0;
%
% Introduce u1 = u_{t} and u2 = u_{x}, then
%       u1_{t} = u2_{x} - u1*u2^2;
%       u2_{t} = u1_{x};
% BCs:  u1(t,0) = 0;    u2(t,1) = 0;
%
% This script tests stability of the PIE, by constructing a distributed
% SOS Lyapunov functional V(u) = <Z(T*u),Pmat*Z(T*u)> + 1e-2*||T*u||^2,
% with Pmat >= 0, and enforcing that the derivative dV of V along
% solutions to the PIE is nonpositive, using piesos_ineq. No ball is
% imposed, so the conditions are enforced globally.
%
% The user may change the degree of the distributed monomials (fourth
% argument of piesos_LF), the options opts of piesos_LF, and the
% coefficient 1e-2 that enforces strict positivity of V.
%
% DJ, 10/15/2026: Initial coding;

clear;  clear stateNameGenerator
pvar x t

% Declare the nonlinear PDE
u1 = pde_var(x,[0,1]);
u2 = pde_var(x,[0,1]);
PDE = [diff(u1,t)==diff(u2,x)-u1*u2^2;
       diff(u2,t)==diff(u1,x);
       subs(u1,x,0)==0;  subs(u2,x,1)==0];

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');  
T = PIE.T;  f = PIE.f;
u = PIE.x;

% Initialize a distributed SOS program
prog = piesos_program(u);

% Declare a distributed SOS candidate Lyapunov functional V, expressed in
% terms of the PDE states T*u, and compute its derivative dV along the PIE
%   V(u) = <Z(T*u),Pmat*Z(T*u)>, Pmat>=0;
%   dV(u) = d/dt V(u), with T*u_t = f(u);
opts.pos = true;       
[prog,V,dV] = piesos_LF(prog,PIE,'pde',1,opts);
% Ensure strict positivity, adding 1e-2*<T*u,T*u>, with derivative
% 2e-2*<T*u,f(u)>
Tu = T*u;      
V = V + 1e-2*innerprod(Tu,Tu); 
dV = dV + 2e-2*innerprod(Tu,f);

% Enforce dV <= -s*g <= 0 along solutions of the PIE
prog = piesos_ineq(prog,-dV);

% Solve the distributed sos program
prog = piesos_solve(prog);
