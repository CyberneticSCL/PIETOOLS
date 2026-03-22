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

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
clear;  clear stateNameGenerator
pvar s t
dom = [0,1];
r = pi^2-0.1;
%r = 1;

% Declare the nonlinear PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+r*x-x*diff(x,s);
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

% Convert ot a PIE
PIE = convert(PDE);
Top = PIE.T;
f = PIE.f;
x = polyopvar(f.varname,s,dom);


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = piesos_program(x);

% Declare monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
Z = dmonomials(x,(1:d)'); % Z_d(v).

% Declare distributed SOS polynomial Vx
% NOTE: Vx is expressed in terms of PDE state
P_opts.deg = 3; % maximal monomial degree of Zop. 
P_opts.exclude = [0;0;0];
P_opts.sep = 1;
[prog,Vx_PDE,Pmat,Zop] = piesos_sosvar(prog,Z,P_opts);
%Vx = quad2lin(Pmat,Zop,Z); % output is in polyopvar

% Ensure strict positivity of the Lyapunov functional.
eppos = 1e-4;
Vx_PDE = Vx_PDE + eppos*innerprod(x,x);
                                                
% Evaluate the derivative of V along the PIE
% NOTEL dV is expressed in terms of PIE state
dV = Liediff(Vx_PDE,PIE); % output is in polyopvar 

% Enforce negativity of the derivative
use_ineq = true;
if use_ineq
    % Use piesos_ineq to enforce negativity
    ineq_opts.psatz = 0:1;
    [prog,W,Qmat,ZQop] = piesos_ineq(prog,-dV,ineq_opts);
else
    % Manually declare distributed SOS polynomial W>=0
    ZQ_degmat = unique(floor(dV.degmat./2),'rows');
    ZQ = polyopvar(f.varname,s,dom);
    ZQ.degmat = ZQ_degmat;
    Q_opts.deg = P_opts.deg+2;
    Q_opts.exclude = [1,0,0]';
    [prog,W,Qmat,ZQop] = piesos_sosvar(prog,ZQ,Q_opts);

    % Enforce dV = -W <= 0
    prog = piesos_eq(prog,dV+W);
end


% Solve the optimization program
prog_sol = piesos_solve(prog);


