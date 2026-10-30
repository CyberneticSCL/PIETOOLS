%
% Burgers eq.
%
% PDE: u_t = u_ss + r*u - u*u_s;
% BCs: u(t,0) = u(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
% maps from fundamental state to PDE states
% u   = (Tv) = \int_0^s T_1(s,t)v(t)dt + \int_s^1 T_2(s,t)v(t)dt;
% u_s = (Rv) = \int_0^s R_1(s,t)v(t)dt + \int_s^1 R_2(s,t)v(t)dt;
% T_1(s,t)=(s-1)*t; T_2(s,t) = s*(t-1);
% R_1(s,t)=t;       R_2(s,t) = t-1;
%
% PIE: (Tv_t) = f(v) = v + (r*Tv) - (Tv)(Rv);
%
% PIE as a polyopvar: Z_1(v) = v; C = C_1 = [[C_111, C_112]
%                                            [C_121, C_122]
%                                            [C_131, C_132]]
%
% C_111 = 1; C_112 = 1; C_121 = 1; C_122 = r*T; C_131 = -T; C_132 = R
% where nz = 1; m_1 = 3; d_1 = d_11 = 2; size(C_1) = (3,2)
%
% This script tests stability of the PIE, by constructing a distributed
% SOS Lyapunov functional V(x) = <Zop*Z(u),Pmat*Zop*Z(u)> + eppos*||u||^2,
% with u = Top*x the PDE state, Z(u) = [u ... u^d]', Zop a basis of
% multiplier and integral operators, and Pmat >= 0. The derivative dV of
% V along solutions to the PIE is enforced to be nonpositive, either using
% piesos_ineq (use_ineq=true), or by declaring a nonnegative distributed
% polynomial functional W and enforcing dV = -W (use_ineq=false). No ball
% is imposed, so the conditions are enforced globally.
%
% The user may change the coefficient r, the degree d, the options
% P_opts (degs, exclude, sep) of the functional, the coefficient eppos,
% and the settings use_ineq, ineq_opts.psatz and Q_opts.degs of the
% constraint on dV.
%
% DJ, 10/15/2026: Initial coding;

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
clear;  clear stateNameGenerator
pvar s t
dom = [0,1];
r = pi^2-0.01;
%r = 1;

% Declare the nonlinear PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+r*x-x*diff(x,s);
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

% Convert ot a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
f = PIE.f;
x = PIE.x;


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = piesos_program(x);

% Declare a distributed SOS polynomial Lyapunov functional V, and compute
% its derivative dV along the PIE.
% V is expressed in terms of the PDE state u = Top*x, with
%   V(x) = <Zop*Z(u),Pmat*Zop*Z(u)>, Pmat>=0, Z(u) = [u ... u^d]';
% dV is the derivative of V along solutions to the PIE, in terms of the
% PIE state x (i.e. Top*x_t = f(x)).
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
P_opts.degs = 3;        % maximal monomial degree of Zop. 
P_opts.exclude = [0;0;0];
P_opts.sep = 1;
P_opts.pos = true;      % V is a distributed SOS polynomial
[prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pde',(1:d)',P_opts);     % all PDE states

% Ensure strict positivity of the Lyapunov functional, adding
% eppos*<u,u>, with derivative 2*eppos*<u,f(x)>.
eppos = 1e-4;
u = Top*x;
Vx = Vx + eppos*innerprod(u,u);
dV = dV + 2*eppos*innerprod(u,f);

% Enforce negativity of the derivative
use_ineq = true;
if use_ineq
    % Use piesos_ineq to enforce negativity
    ineq_opts.psatz = 1;    % multipliers (b-s), (s-a), see piesos_quadopvar
    [prog,W,Qmat,ZQop] = piesos_ineq(prog,-dV,ineq_opts);
else
    % Manually declare distributed SOS polynomial W>=0
    % dV has unused variables besides the PIE state x, 
    % keep only the columns of x
    [~,idcs_x] = ismember(x.varname,dV.varname); 
    ZQ_degmat = unique(floor(dV.degmat(:,idcs_x)./2),'rows');
    ZQ = dmonomials(x,ZQ_degmat);   % degmat is not settable
    Q_opts.degs = P_opts.degs+2;
    Q_opts.exclude = [1,0,0]';
    [prog,W,Qmat,ZQop] = piesos_sosvar(prog,ZQ,Q_opts);

    % Enforce dV = -W <= 0
    prog = piesos_eq(prog,dV+W);
end

% Solve the optimization program
prog_sol = piesos_solve(prog);
