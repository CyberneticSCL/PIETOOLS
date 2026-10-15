%
% Coupled Burgers equations
%
% PDE: u1_t = u1_ss + r*u1 - u1*u1_s;
%      u2_t = u2_ss + r*u2 - u2*u2_s;
% BCs: u1(t,0) = 0;     u1(t,0.5) = u2(t,0);
%      u2(t,0.5) = 0;   u1_{s}(t,0.5) = u2_{s}(t,0);
%
% fundamental state: v = [v1;v2] = [u1_ss; u2_ss] \in L_2[0,0.5];
%
% PIE: (T(1,:)*v_t) = f1(v) = v1 + (r*T(1,:)*v) - (T(1,:)*v)(R(1,:)*v);
%      (T(2,:)*v_t) = f2(v) = v2 + (r*T(2,:)*v) - (T(2,:)*v)(R(2,:)*v);
% where u = T*v and u_s = R*v for u = [u1;u2].
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
% Popts (degs, exclude, sep) of the functional, the coefficient eppos,
% and the settings use_ineq, ineq_opts.psatz, Q_opts.degs and
% Q_opts.psatz of the constraint on dV.
%
% DJ, 10/15/2026: Initial coding;

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
clear;  clear stateNameGenerator
pvar s s_dum t
dom = [0,0.5];
r = pi^2-0.01;

% Declare the nonlinear PDE
x1 = pde_var(s,dom);
x2 = pde_var(s,dom);
PDE = [diff(x1,t)==diff(x1,s,2)+r*x1-x1*diff(x1,s);
       diff(x2,t)==diff(x2,s,2)+r*x2-x2*diff(x2,s);
       subs(x1,s,dom(1))==0;     subs(x1,s,dom(2))==subs(x2,s,dom(1));
       subs(diff(x1,s),s,dom(2))==subs(diff(x2,s),s,dom(1));
       subs(x2,s,dom(2))==0];

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
f = PIE.f;

%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
vartab = PIE.x;
prog = piesos_program(vartab);

% Declare a distributed SOS Lyapunov functional V, and compute its
% derivative dV along the PIE.
% V is expressed in terms of the PDE states u = Top*x, with
%   V(x) = <Zop*Z(u),Pmat*Zop*Z(u)>, Pmat>=0, Z(u) = [u ... u^d]';
% dV is the derivative of V along solutions to the PIE, in terms of the
% PIE state x (i.e. Top*x_t = f(x)).
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
Popts.degs = 2;         % maximal monomial degree of Zop.
Popts.exclude = [0;0;0];
Popts.sep = true;
Popts.pos = true;       % V is a distributed SOS polynomial
[prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pde',(1:d)',Popts);

% Ensure strict positivity of the Lyapunov functional, adding
% eppos*<u,u>, with derivative 2*eppos*<u,f(x)>.
eppos = 1e-4;
u = Top*vartab;
Vx = Vx + eppos*innerprod(u,u);
dV = dV + 2*eppos*innerprod(u,f);

use_ineq = false;
if use_ineq
    % Use piesos_ineq to enforce
    ineq_opts.psatz = 1;    % multipliers (b-s), (s-a), see piesos_quadopvar
    [prog,W,Qmat,ZQop] = piesos_ineq(prog,-dV,ineq_opts);
else
    % Declare a nonnegative distributed polynomial functional W
    Q_opts.degs = Popts.degs+4;
    % dV has unused variables besides the PIE states x, keep only the columns of x
    [~,idcs_x] = ismember(vartab.varname,dV.varname);
    ZQ_degmat = unique(floor(dV.degmat(:,idcs_x)./2),'rows');
    ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
    ZQ = dmonomials(vartab,ZQ_degmat);
    Q_opts.exclude = [1,0,0]';
    Q_opts.psatz = 1;       % multipliers (b-s), (s-a), see piesos_quadopvar
    [prog,W,Qmat,ZQop] = piesos_sosvar(prog,ZQ,Q_opts);

    % Enforce dV = -W <= 0
    prog = piesos_eq(prog,dV+W);
end

% Solve the optimization program
prog_sol = piesos_solve(prog);
