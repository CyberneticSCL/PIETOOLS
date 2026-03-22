%
% Kuramoto-Sivashinsky Equation
%
% PDE: u_t = -u_ssss - u_ss - -r*u*u -u*u_s;
% BCs: u(t,0) = 0;     u1(t,1) = 0;
%      u_s(t,0) = 0;   u1_{s}(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
%
% PIE: (T*v_t) = f1(v) = -v -R2*v - r*(T*v)*(T*v) -(T*v)*(R*v);
%

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
pvar s s_dum t
dom = [0,1];
r = 0.5;           % should be [−0.6500, 0.7202]

% Declare a PDE and convert to PIE to get the relevant maps from PIE state
% to PDE state
x = pde_var(s,dom);
PDE = [diff(x,t)==-diff(x,s,4)-diff(x,s,2)-r*x^2-x*diff(x,s);
        subs(x,s,0)==0;     subs(diff(x,s),s,0)==0;
        subs(x,s,1)==0;     subs(diff(x,s),s,1)==0];

% Convert to PIE
PIE = convert(PDE);
Top = PIE.T;
f = PIE.f;




%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
vartab = polyopvar(f.varname,s,dom); % Z_d(v).
prog = piesos_program(vartab);

% Declare monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
Z = dmonomials(vartab,(1:d));

% Declare distributed SOS polynomial Vx
% NOTE: Vx is expressed in terms of PDE state
pdegs = 7; % maximal monomial degree of Zop. 
Popts.deg = pdegs;
Popts.exclude = [0;1;1];
Popts.sep = true;
[prog,Vx_PDE,Pmat,Zop] = piesos_sosvar(prog,Z,Popts);
% [prog,Pmat,Zop] = piesos_poslpivar(prog,Z,pdegs,Popts);
% Vx = quad2lin(Pmat,Zop,Z);          

% Ensure strict positivity of the Lyapunov functional.
eppos = 1e-4;
Vx_PDE = Vx_PDE + eppos*innerprod(vartab,vartab);
                                                
% Evaluate the derivative of V along the PIE
% NOTEL dV is expressed in terms of PIE state
dV = Liediff(Vx_PDE,PIE); % output is in polyopvar 

% Enforce negativity of the derivative
use_ineq = true;
if use_ineq
    ineq_opts.psatz = 0;
    [prog,W,Qmat,ZQop] = piesos_ineq(prog,-dV,ineq_opts);
else
    % Declare a nonnegative distributed polynomial functional W
    qdegs = pdegs+4;
    ZQ_degmat = unique(floor(dV.degmat./2),'rows');
    ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
    ZQ = vartab;
    ZQ.degmat = ZQ_degmat;
    Q_opts.deg = qdegs;
    Q_opts.exclude = [1,0,0]';
    Q_opts.psatz = 0;
    [prog,W,Qmat,ZQop] = piesos_sosvar(prog,ZQ,Q_opts);
    % [prog,Qmat,ZQop] = piesos_poslpivar(prog,ZQ,qdegs,Q_opts);
    % W = quad2lin(Qmat,ZQop,ZQ);

    % Enforce dV = -W <= 0
    prog = piesos_eq(prog,dV+W);
end


% Solve the optimization program
solve_opts.simplify = true;
prog_sol = piesos_solve(prog);


