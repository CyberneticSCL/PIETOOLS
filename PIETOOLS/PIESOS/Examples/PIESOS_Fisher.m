%
% Fisher eq.
%
% PDE: u_t = u_ss + alp*u - bet*u^2;
% BCs: u(t,0) = u(t,1) = 0;
%
% fundamental state: v = u_ss \in L_2[0,1];
%
% PIE: (Tv_t) = f(v) = v + alp*(Tv) - bet*(Tv)^2;
%
% For alp=5 and bet=-1, the PDE is locally stable on a ball in the L2
% norm of radius R ~ 4.405.
%
% This script tests local stability of the PIE on the L2 ball
% g = R^2 - ||u||^2 >= 0, by constructing a distributed SOS Lyapunov
% functional V(u) = <Z(u),P*Z(u)> + eppos*||u||^2, with u = T*v the PDE
% state, Z(u) = [u ... u^d]', and P a multiplier operator, with
% P(s) = Z(s)'*Q*Z(s) and Q >= 0. It enforces dV <= 0 on the set
% G = {u | g(u) >= 0}, using an SOS multiplier lam2 of degree d_psatz2
% (or globally, if R <= 0), and, if use_bnd=true, the upper bound
% V <= gam*||u||^2, minimizing gam. The positivity constraints are
% enforced with piesos_ineq if use_ineq1 or use_ineq2 is true, and by
% explicitly declaring distributed SOS functionals otherwise.
%
% The user may change the coefficients alp and bet (section 1), and the
% settings in section 2: the degrees d, d_psatz1, d_psatz2 and pdeg, the
% coefficient eppos, the decay rate k, the
% radius R, and the flags use_ineq1, use_ineq2 and use_bnd.
%
% DJ, 10/15/2026: Initial coding;

%%%% 1. Declare the PIE for which to test stability

% Declare spatial variables, domain, and parameters.
clear;  clear stateNameGenerator
pvar s t
dom = [0,1];
alp = 5;
bet = -1;

% Declare the nonlinear PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+alp*x-bet*x^2;
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
f = PIE.f;
x = PIE.x;


%%%% 2. Test stability using PIESOS

% % Declare settings
% Set distributed monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
% Set distributed monomial basis for psatz multiplier
d_psatz1 = 0;
Zg1 = dmonomials(x,(1:d_psatz1));
d_psatz2 = 1;
Zg2 = dmonomials(x,(1:d_psatz2));
% Monomial degree in independent variables, used to parameterize Pop
pdeg = 4;
% Coercivity of Pop
eppos = 1;
% Enforce positivity using piesos_ineq
use_ineq1 = false;
use_ineq2 = false;
% Enforce strict negativity of LF derivative
k = 0;              % Rate of exponential decay
% Enforce local decay on semialgebraic set where g>=0
R = 4.0;         % feasible up to R~4.0479
g = R^2-innerprod(Top*x,Top*x);                             % sphere in L2 norm
% Enforce an upper bound on the LF as V(x) <= M*||T*x||_{L2}^2
use_bnd = true;


% % Initialize PIESOS program structure.
if use_bnd
    dpvar gam
    prog = piesos_program(x,gam);
    % Minimize the upper bound gamma on the Lyapunov function
    prog = piesos_setobj(prog,gam);
    %prog = piesos_ineq(prog,gam);
else
    prog = piesos_program(x);
end


% % Declare the Lyapunov functional and its derivative.
%       V(u) = <u, P u>_{L2}, u = T*x,
%   where P is a multiplier operator, P(s) = Z(s)'*Q*Z(s) >= 0. Then,
%       d/dt <u, P u>_{L2} = 2*<T*x, P f(x)>_{L2}

% % First, parameterize an SOS multiplier operator, and compute the
% % derivative of the functional.
disp(" --- Declaring the distributed SOS Lyapunov functional ---")

% Declare V and its derivative dV in terms of the PIE state x.
Popts.degs = pdeg;   
Popts.exclude = [0;1;1];    % multiplier operator only 
Popts.pos = true; 
[prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pde',(1:d)',Popts);
% Ensure strict positivity of the Lyapunov functional, adding
% eppos*<u,u>, with derivative 2*eppos*<u,f(x)>.
Tx = Top*x;
Vx = Vx + eppos*innerprod(Tx,Tx); 
dV = dV + 2*eppos*innerprod(Tx,f);


% % Add the psatz term, to allow dV>0 when g = R-<Tx,Tx> <=0
if R>0
    disp(" --- Declaring the psatz term ---")
    % Declare the multiplier for the upper bound
    if use_bnd
        if d_psatz1>=1
            lam1_opts.exclude = [1,0,0]';
            lam1_opts.degs = pdeg;    
            lam1_opts.psatz = 0;
            [prog,lam1] = piesos_sosvar(prog,Zg1,lam1_opts);
            V_bnd = gam*innerprod(Tx,Tx)-Vx-lam1*g;  %(gam*||Tx||^2 - Vx) >= lam1*g >= 0
        else
            V_bnd = gam*innerprod(Tx,Tx)-Vx;
        end
    end

    % Declare the multiplier for the derivative
    lam2_opts.exclude = [1,0,0]';
    lam2_opts.degs = pdeg;
    lam2_opts.psatz = 0;
    [prog,lam2] = piesos_sosvar(prog,Zg2,lam2_opts);
    Wgg = lam2*g;
    dV_g = dV+Wgg;
else
    % Global stability test
    dV_g = dV;
end


% % Enforce Mx-Vx >= 0 using distributed SOS
if use_bnd
    disp(" --- Enforcing the bound on the Lyapunov functional ---")
    if use_ineq1
        disp("  --  using piesos_ineq")
        ineq1_opts.psatz = 0;
        [prog,W1,Q1mat,ZQ1op] = piesos_ineq(prog,V_bnd,ineq1_opts);
    else
        % Declare W >= 0
        disp("  --  manually declaring a distributed SOS functional")
        % Determine monomials to balance those of dV
        % V_bnd has unused variables besides the PIE state x, keep only the columns of x
        [~,idcs_x] = ismember(x.varname,V_bnd.varname);    
        Z_bnd_degmat = unique(floor(V_bnd.degmat(:,idcs_x)./2),'rows'); 
        Z_bnd = dmonomials(x,Z_bnd_degmat);                         
        % Declare distributed SOS functional in terms of these monomials
        Q1_opts.degs = pdeg;          
        Q1_opts.exclude = [1,0,0]';
        Q1_opts.psatz = 1;      % multipliers (b-s), (s-a)      
        [prog,W1,Q1mat,ZQ1op] = piesos_sosvar(prog,Z_bnd,Q1_opts);

        % Enforce gam*||u||^2-V = W >= 0
        disp("  --  enforcing equality")
        prog = piesos_eq(prog,V_bnd-W1);
    end
end


% % Enforce dV <= 0 using distributed SOS
disp(" --- Enforcing negativity of the derivative ---")
if use_ineq2
    disp("  --  using piesos_ineq")
    ineq2_opts.psatz = 0;
    [prog,W2,Q2mat,ZQ2op] = piesos_ineq(prog,-dV_g-2*k*Vx,ineq2_opts);
else
    % Declare W >= 0
    disp("  --  manually declaring a distributed SOS functional")
    % Determine monomials to balance those of dV
    % dV_g has unused variables besides the PIE state x, keep only the columns of x
    [~,idcs_x] = ismember(x.varname,dV_g.varname);    
    ZQ_degmat = unique(floor(dV_g.degmat(:,idcs_x)./2),'rows'); 
    ZQ = dmonomials(x,ZQ_degmat); 
    % Declare distributed SOS functional in terms of these monomials
    qdeg = pdeg+1;
    if d_psatz2==2
        Q2_opts.degs = {[qdeg,qdeg]; [qdeg,round(qdeg/2)]; [qdeg-2,round(qdeg/4)]};
    else
        Q2_opts.degs = qdeg; 
    end
    % Q2_opts.deg = pdeg+1;
    Q2_opts.exclude = [1,0,0]';
    Q2_opts.psatz = 1;      % multipliers (b-s), (s-a)      
    [prog,W2,Q2mat,ZQ2op] = piesos_sosvar(prog,ZQ,Q2_opts);

    % Enforce dV = -W <= 0
    disp("  --  enforcing equality")
    prog = piesos_eq(prog,dV_g+2*k*Vx+W2);
end


% % Solve the optimization program
disp(" --- Solving the optimization program ---")
sol_opts.simplify = true;                           
prog_sol = piesos_solve(prog,sol_opts); 


% % Extract the solution
sol_info = prog_sol.solinfo.info;
if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
    disp("PIESOS program was not solved.")
else
    disp("PIESOS program was successfully solved.")
    Vsol = piesos_getsol(prog_sol,Vx);
    if use_bnd
        gam_sol = piesos_getsol(prog_sol,gam);
    end
end
