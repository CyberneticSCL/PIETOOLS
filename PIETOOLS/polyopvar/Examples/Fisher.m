%
% Fishers eq.
%
% PDE: u_t = u_ss + alp*u - bet*u^2;
% BCs: u(t,0) = u(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
%
% PIE: (Tv_t) = f(v) = v + (alp*Tv) - bet*(Tv)^2;
%
%

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
clear;  clear stateNameGenerator
pvar s t
dom = [0,1];
alp = -0.1;
bet = 0.15;
R = 1;          % Radius of ball in which to test stability

% Declare the nonlinear PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+alp*x-bet*x^2;
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

% Convert ot a PIE
PIE = convert(PDE);
Top = PIE.T;
f = PIE.f;


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = lpiprogram(s,dom);

% Declare monomial basis of SOS Lyapunov functional.
x = polyopvar(f.varname,s,dom);
d = 2;                  % degree of distributed monomial basis, will be doubled in LF
Z = dmonomials(x,(1:d)); % Z_d(v).

% Declare PSD operator acting on degree d monomial basis and add variables to LPI program.
% opts = ;
P_opts.deg = 0; % maximal monomial degree of Zop. 
P_opts.exclude = [0;0;0];
P_opts.sep = 1;
[prog,Vx,Pmat,Zop] = piesos_sosvar(prog,Z,P_opts);

% Ensure strict positivity of the Lyapunov functional.
eppos = 1e-4;
Vx = Vx + eppos*innerprod(x,x);
                                                
% Evaluate the derivative of V along the PIE
dV = Liediff(Vx,PIE); % output is in polyopvar 

% Declare a nonnegative distributed polynomial functional W
ZQ_degmat = unique(floor(dV.degmat./2),'rows');
ZQ = polyopvar(f.varname,s,dom);
ZQ.degmat = ZQ_degmat;
Q_opts.deg = P_opts.deg+2;
Q_opts.exclude = [1,0,0]';
[prog,W,Qmat,ZQop] = piesos_sosvar(prog,ZQ,Q_opts);

% Allow W<=0 when g = R-<Tx,Tx> <=0
g = R-innerprod(Top*x,Top*x);
Wg_opts = Q_opts;
Zg_degmat = unique(floor((W.degmat-max(g.degmat,[],1))./2),'rows');
Zg_degmat = Zg_degmat(~any(Zg_degmat<0,2),:);
Zg = ZQ;
if any(Zg_degmat>0)
    Zg_degmat = Zg_degmat(~all(Zg_degmat,2)==0,:);
end
Zg.degmat = Zg_degmat;
[prog,Wg] = piesos_sosvar(prog,Zg,Wg_opts);
W = W + Wg*g;

% Enforce dV = -W <= 0
prog = piesos_eq(prog,dV+W);

% Solve the optimization program
prog_sol = lpisolve(prog);

