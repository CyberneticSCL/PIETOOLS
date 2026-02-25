%
% Coupled Burgers equations
%
% PDE: u1_t = u1_ss + r*u1 - u1*u1_s;
%      u2_t = u2_ss + r*u2 - u2*u2_s;
% BCs: u1(t,0) = 0;     u1(t,0.5) = u2(t,0);
%      u2(t,0.5) = 0;   u1_{s}(t,0.5) = u2_{s}(t,0);
%
% fundamental state: v = [v1;v2] = [u1_ss; u2_ss] \in L_2[0,1];
%
%
% PIE: (T(1,:)*v_t) = f1(v) = v1 + (r*T(1,:)*v) - (T(1,:)*v)(R(1,:)*v);
%      (T(2,:)*v_t) = f2(v) = v2 + (r*T(2,:)*v) - (T(2,:)*v)(R(2,:)*v);
%

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
clear;  clear stateNameGenerator
pvar s s_dum t
dom = [0,0.5];
r = pi^2-0.1;
%r = 1;

% Declare the nonlinear PDE
x1 = pde_var(s,dom);
x2 = pde_var(s,dom);
PDE = [diff(x1,t)==diff(x1,s,2)+r*x1-x1*diff(x1,s);
       diff(x2,t)==diff(x2,s,2)+r*x2-x2*diff(x2,s);
       subs(x1,s,dom(1))==0;     subs(x1,s,dom(2))==subs(x2,s,dom(1));
       subs(diff(x1,s),s,dom(2))==subs(diff(x2,s),s,dom(1));
       subs(x2,s,dom(2))==0];

% Convert ot a PIE
PIE = convert(PDE);
T1op = PIE.T;
f = PIE.f;


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = lpiprogram([s,s_dum],dom);

% Declare monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
vartab = polyopvar(f.varname,s,dom); % Z_d(v).
Z = dmonomials(vartab,(1:d));

% Declare PSD operator acting on degree d monomial basis and add variables to LPI program.
% opts = ;
pdegs = 3; % maximal monomial degree of Zop. 
Popts.exclude = [0;0;0];
Popts.sep = true;
[prog,Pmat,Zop] = soslpivar(prog,Z,pdegs,Popts);
Vx = quad2lin(Pmat,Zop,Z);          

% Add PD constraint to LPI program.
eppos = 1e-4;
is_diag_term = sum(Vx.degmat,2)==2 & max(Vx.degmat,[],2)==2;
for i=find(is_diag_term')
    idx = find(Vx.C.ops{i}.omat(1,1)==0);
    Vx.C.ops{i}.params(idx) = Vx.C.ops{i}.params(idx) + eppos;  % can be done more elegantly once polyopvar is better developed
end

% Evaluate the derivative of V along the PIE
dV = Liediff(Vx,PIE);

% Declare a nonnegative distributed polynomial functional W
qdegs = pdegs+2;
ZQ_degmat = unique(floor(dV.degmat./2),'rows');
ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
ZQ = vartab;
ZQ.degmat = ZQ_degmat;
Q_opts.exclude = [1,0,0]';
Q_opts.psatz = 0;
[prog,Qmat,ZQop] = soslpivar(prog,ZQ,qdegs,Q_opts);
W = quad2lin(Qmat,ZQop,ZQ);

% Enforce dV = -W <= 0
prog = soslpi_eq(prog,dV+W);

% Solve the optimization program
prog_sol = lpisolve(prog);
