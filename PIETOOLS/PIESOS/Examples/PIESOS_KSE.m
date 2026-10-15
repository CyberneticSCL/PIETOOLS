%
% Kuramoto-Sivashinsky Equation
%
% PDE: u_t = -u_ssss - u_ss - r*u*u - u*u_s;
% BCs: u(t,0) = 0;     u(t,1) = 0;
%      u_s(t,0) = 0;   u_s(t,1) = 0;
%
% fundamental state: v = u_ssss \in L_2[0,1];
%
% PIE: (T*v_t) = f(v) = -v - R2*v - r*(T*v)*(T*v) - (T*v)*(R*v);
% where u = T*v, u_s = R*v and u_ss = R2*v.
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
% The user may change the coefficient r, the degree d, the maximal
% monomial degree pdegs and options Popts (exclude, sep) of the
% functional, the coefficient eppos, and the settings use_ineq,
% ineq_opts.psatz, qdegs and Q_opts.psatz of the constraint on dV.
%
% DJ, 10/15/2026: Initial coding;

%%%% 1. Declare the PIE for which to verify stability

% Declare spatial variables, domain, and parameters.
pvar s s_dum t
dom = [0,1];
r = 0.3;           % should be [-0.6500, 0.7202]

% Declare the PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==-diff(x,s,4)-diff(x,s,2)-r*x^2-x*diff(x,s);
        subs(x,s,0)==0;     subs(diff(x,s),s,0)==0;
        subs(x,s,1)==0;     subs(diff(x,s),s,1)==0];

% Convert to PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
f = PIE.f;


%%%% 2. Solve distributed SOS program for stability

% Set up PIESOS program structure.
vartab = PIE.x;
prog = piesos_program(vartab);

% Declare a distributed SOS Lyapunov functional V, and compute its
% derivative dV along the PIE.
% V is expressed in terms of the PDE state u = Top*x, with
%   V(x) = <Zop*Z(u),Pmat*Zop*Z(u)>, Pmat>=0, Z(u) = [u ... u^d]';
% dV is the derivative of V along solutions to the PIE, in terms of the
% PIE state x (i.e. Top*x_t = f(x)).
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
pdegs = 7; % maximal monomial degree of Zop.
Popts.degs = pdegs; 
Popts.exclude = [0;1;1];
Popts.sep = true;
Popts.pos = true;       % V is a distributed SOS polynomial
[prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pde',(1:d)',Popts);

% Ensure strict positivity of the Lyapunov functional, adding
% eppos*<u,u>, with derivative 2*eppos*<u,f(x)>.
eppos = 1e-2;
u = Top*vartab;
Vx = Vx + eppos*innerprod(u,u);
dV = dV + 2*eppos*innerprod(u,f);

% Enforce negativity of the derivative
use_ineq = false;
if use_ineq
    ineq_opts.psatz = 0;
    [prog,W,Qmat,ZQop] = piesos_ineq(prog,-dV,ineq_opts);
else
    % Declare a nonnegative distributed polynomial functional W
    qdegs = pdegs+5;
    % dV has unused variables besides the PIE state x, keep only the columns of x
    [~,idcs_x] = ismember(vartab.varname,dV.varname);
    ZQ_degmat = unique(floor(dV.degmat(:,idcs_x)./2),'rows');
    ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
    ZQ = dmonomials(vartab,ZQ_degmat); 
    Q_opts.degs = qdegs; 
    Q_opts.exclude = [1,0,0]';
    Q_opts.psatz = 0;
    [prog,W,Qmat,ZQop] = piesos_sosvar(prog,ZQ,Q_opts);

    % Enforce dV = -W <= 0
    prog = piesos_eq(prog,dV+W);
end

% Solve the optimization program
prog_sol = piesos_solve(prog);
