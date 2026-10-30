%
% Heat equation with nonlinearity and stabilizing feedback
%
% PDE: u_t = a*u_ss + u^2 - k/L * \int_0^L u ds;
% BCs: u_s(t,0) = u(t,L) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1] (assuming L=1);
%
% PIE: (Tv_t) = f(v) = (a - k/L * vec(T))*v + Tv*Tv;
%
% This script verifies decay of the functional
%   V(u) = <u_{s},u_{s}> = -<x,T*x>
% on the region
%   G:={u in D| <u_{s},u_{s}><=q^2},
% where x is the PIE state, by maximizing the decay rate epneg in
%   dV <= -epneg*||x||^2   on G,
% with dV the derivative of V along solutions to the PIE. This constraint
% is enforced using piesos_ineq, with an SOS multiplier lam1 of degree
% d_psatz for the region.
%
% The user may change the parameters in "Script parameters" and the
% coefficients L, a and k of the PDE: the radius q of G, the degree
% d_psatz of the multiplier, the degree lam_deg of its monomials, and the
% options ineq_opts (psatz, degs) of piesos_ineq.
%
% DJ, 10/15/2026: Initial coding;
% DJ, 10/29/2026: Use the regulated output polynomial PIE.gz, instead of PIE.h;

clear;  clear stateNameGenerator

% Declare the nonlinear PDE
pvar   s t s_dum
L = 1;
a = 1;
k = 2*a*pi^2/L^2;
dom = [0,L];
x = pde_var(s,dom);

% FULL PDE
z = pde_var('out',s,dom);    % output z = u_s, to obtain Rop: z = Rx  
PDE = [diff(x,t) == a*diff(x,s,2) +x^2- (k/L)*int(x,s,dom(1),dom(2));    
       subs(diff(x,s,1),s,dom(1))==0; subs(x,s,dom(2))==0;
       z==diff(x,s)];                                             

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
ux = PIE.gz;     % u_{x}
f = PIE.f;
x = PIE.x;  

% Script parameters
q = 1;              % for k=2*a*pi^2/L^2, the set ||u_{x}||_{L2}<= a0*pi^2/(2*L^(5/2)) should be an invariant region

d_psatz = 1;    % distributed monomial degree of multiplier for regional stability

lam_deg = 2;    % degree of kernel monomials used to define the operator lam1
ineq_opts.psatz = 1;
ineq_opts.degs = 4;     % degree of kernel monomials in negativity constraint

%% Declare the stability PIESOS program

% Declare the LF
%   V(x) = -<x,Tx> = <u_{x},u_{x}>
V = -innerprod(x,Top*x);

% Compute the Lie derivative
dV = -2*innerprod(x,f);

% Enforce negativity on the specified region
g = q^2 - innerprod(ux,ux);

% Maximize the decay rate
dpvar epneg
prog = piesos_program(x,epneg);
prog = piesos_setobj(prog,-epneg);
%prog = piesos_ineq(prog,epneg);
Zg1 = dmonomials(x,(1:d_psatz));
lam1_opts.exclude = [1,0,0]';
lam1_opts.degs = lam_deg;          % degree of monomials (not distributed) used to define the operator lam1
lam1_opts.psatz = 0;
[prog,lam1] = piesos_sosvar(prog,Zg1,lam1_opts);

prog = piesos_ineq(prog,-dV-lam1*g-epneg*innerprod(x,x),ineq_opts);

% Solve the program
sol_opts.simplify = true;
prog_sol = piesos_solve(prog,sol_opts);
eps_opt = double(piesos_getsol(prog_sol,epneg));
if eps_opt<0
    fprintf(2,"\n Warning: stability was verified with negative decay rate.\n")
end

