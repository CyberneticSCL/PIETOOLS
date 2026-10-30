%
% Fisher eq.
%
% PDE: u_t = u_ss + alp*u - bet*u^2;
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
% PIE: (Tv_t) = f(v) = v + alp*(Tv) - bet*(Tv)^2;
%
% This script tests local exponential stability of the PIE, by
% constructing a distributed SOS Lyapunov functional
% V(v) = <Zop*Z(v),Pmat*Zop*Z(Tv)>, with Z(v) = [v ... v^d]' and Zop a
% basis of multiplier and integral operators, such that
%   V >= eppos*||u||^2;         (lower bound)
%   V <= gam*||u||^2;           (upper bound, gam is minimized)
%   dV <= -2*k*V;               (decay at rate k along the PIE)
% The upper bound and the decay condition are enforced on the ball
% ||u||^2 <= R^2 if the degrees d_psatz1 and d_psatz2 of the SOS
% multipliers are at least 1, and globally if they are 0.
% If BALL=false, the norm ||u|| in the bounds is replaced by ||u_s||, and
% the ball is the Sobolev ball ||u||^2 + ||u_s||^2 <= R^2.
%
% The user may change the parameters in "Script parameters": the
% coefficients alp and bet, the radius R, the decay rate k, the degrees d
% and pdeg of the functional, BALL, and the degrees d_psatz1, d_psatz2,
% q1_deg, q2_deg, q3_deg, lam1_deg and lam2_deg of the SOS multipliers
% and variables.
%
% Remark: for the Burgers-Fisher eq. u_t = u_ss + alp*u - u*u_s - bet*u^2,
% the zero solution should be globally stable when
%   (1/pi^2)*exp(-3*bet)*(alp+(9/2)*bet^2) < 1,
% or perhaps, using a weighted Poincare inequality, when
%   alp+(9/4)*bet^2 < 1,
% using V(u) = \int_0^1 exp(-3*bet*s)*u(s)^2 ds.
%
% DJ, 10/15/2026: Initial coding;
% DJ, 10/29/2026: Use the regulated output polynomial PIE.gz, instead of PIE.h;
%

clear;  clear stateNameGenerator

% Declare the nonlinear PDE
pvar   s t s_dum
dom =  [0,1];
x   = pde_var(s,dom);
alp =  5;
bet = -1;
z = pde_var('out',s,dom);    % output z = u_s, to obtain Rop: z = Rx
PDE = [diff(x,t)==diff(x,s,2)+alp*x-bet*x^2;
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0;
       z==diff(x,s)];


% Script parameters
R = 2.95; % Radius of ball in which to test stability - feasible up to R~4.0479.
k = 1.3;    % Rate of decay of the functional.
d = 1;    % degree of LF distributed monomial basis (will be doubled in LF).
pdeg = 4; % degree of Zs monomials of positive P operator.
BALL = true; % local test on L2 ball (if TRUE) or Sobolev ball (if false) of radius R.
d_psatz1 = 1; % degree of distributed monomial in upper bound condition of LF.
d_psatz2 = 1; % degree of distributed monomial in upper bound condition of LF derivative.
eppos = 1.0; % coefficient in LF lower bound condition.
q1_deg = 4; % degree of monomials (not distributed) used to define the operator W1
q2_deg = 4; % degree of monomials (not distributed) used to define the operator W2
q3_deg = 4; % degree of monomials (not distributed) used to define the operator W3
lam1_deg = 4; % degree of monomials (not distributed) used to define the operator lam1
lam2_deg = 4; % degree of monomials (not distributed) used to define the operator lam2

%% 1. Modelling PIE.

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
% The Sobolev ball requires Rx = Rop*x, the derivative u_s in terms of the PIE state
% x. We obtain it by declaring an output z = u_s in the PDE, which the converter
% returns as the 'polyopvar' PIE.gz = Rop*x:
f = PIE.f;
x = PIE.x;


%% 2. Setting up PIESOS program

% Initialize PIESOS program structure.
dpvar gam
prog = piesos_program(x,gam);

% Define as a minimization problem on the upper bound (gamma) of the LF
% over the domain L_{2,R}.
prog = piesos_setobj(prog,gam);

%% 3. Construct LF
% V(v) = <Z_d(v), Pop*Z_d(Tv)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Declare the Lyapunov functional.
disp(" --- Declaring the distributed SOS Lyapunov functional ---")

% Declare the functional, and its derivative along the PIE, in terms of the
% PIE state, with
%   V(x) = <Zop*Z(x), Pmat*Zop*Z(Tx)>,
% for Pmat constant (not enforcing positivity here), and Zop a basis of
% multiplier and integral operators with monomials of degree up to pdeg.
% The symmetry condition V(xL,T*xR) = V(xR,T*xL) is enforced by piesos_LF.
Popts.degs = struct('in',pdeg,'out',0);
Popts.exclude = [0;0;0];
[prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pie',(1:d)',Popts);
Tx = Top*x;

%% 4. Construct the derivative of the LF along the PIE (for up to degree 3).
% The derivative dV was computed together with V above.
disp(" --- Computing the derivative along the PIE ---")


%% 5. Define the specified local ball as a semialgebraic set.
disp(" --- Declaring the local ball ---")
if BALL
    g = R^2-innerprod(Tx,Tx); % L2 ball of radius R.
else
    Rx = PIE.gz;
    g = R^2-(innerprod(Top*x,Top*x) + innerprod(Rx,Rx)); % Sobolev ball of radius R.
end

%% 6. Define the lower bound on the LF (holds globally) and enforce constraint.

if BALL % LF upper bound dependent on ball.
    V_low = Vx - eppos*innerprod(Tx,Tx);
else
    V_low = Vx - eppos*innerprod(Rx,Rx);
end

% Enforce lower bound by defining a SOS distributed monomial
disp(" --- Enforcing the bound on the Lyapunov functional ---")

% Declare W1 >= 0
% V_low has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,V_low.varname);
Z_bnd_degmat = unique(floor(V_low.degmat(:,idcs_x)./2),'rows'); % This is the degree of W1
Z_bnd = dmonomials(x,Z_bnd_degmat);
Q1_opts.degs = q1_deg; % degree of monomials (not distributed) used to define the operator W1
Q1_opts.exclude = [1,0,0]';
Q1_opts.psatz = 1;      % multipliers (b-s), (s-a)
[prog,W1,Q1mat,ZQ1op] = piesos_sosvar(prog,Z_bnd,Q1_opts);

disp("  --  enforcing lower bound equality")
prog = piesos_eq(prog,V_low-W1);

%% 7. Define the upper bound on the LF over the specified ball and enforce constraint.

% Declare the SOS multiplier, lam1, then define bound.
Zg1 = dmonomials(x,(1:d_psatz1));
if d_psatz1>=1
    lam1_opts.exclude = [1,0,0]';
    lam1_opts.degs = lam1_deg; % degree of monomials (not distributed) used to define the operator lam1
    lam1_opts.psatz = 0;
    [prog,lam1] = piesos_sosvar(prog,Zg1,lam1_opts);

    if BALL % LF upper bound dependent on ball.
        V_up = gam*innerprod(Tx,Tx)-Vx-lam1*g;
    else
        V_up = gam*innerprod(Rx,Rx)-Vx-lam1*g;
    end
else
    if BALL % LF upper bound dependent on ball.
        V_up = gam*innerprod(Tx,Tx)-Vx;
    else
        V_up = gam*innerprod(Rx,Rx)-Vx;
    end
end

% Enforce upper bound by defining a SOS distributed monomial
disp(" --- Enforcing the bound on the Lyapunov functional ---")

% Declare W2 >= 0
% V_up has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,V_up.varname);
Z_bnd_degmat = unique(floor(V_up.degmat(:,idcs_x)./2),'rows');
Z_bnd = dmonomials(x,Z_bnd_degmat);
Q2_opts.degs = q2_deg; % degree of monomials (not distributed) used to define the operator W2
Q2_opts.exclude = [1,0,0]';
Q2_opts.psatz = 1;      % multipliers (b-s), (s-a)
[prog,W2,Q2mat,ZQ2op] = piesos_sosvar(prog,Z_bnd,Q2_opts);

disp("  --  enforcing upper bound equality")
prog = piesos_eq(prog,V_up-W2);

%% 8. Enforce symmetry condition
% The symmetry condition V(xL,T*xR) = V(xR,T*xL) is enforced by piesos_LF.

%% 9. Define the upper bound on the LF derivative over the specified ball and enforce constraint.

% Declare the SOS multiplier, lam2, then define bound.
Zg2 = dmonomials(x,(1:d_psatz2));
if d_psatz2>=1
    lam2_opts.exclude = [1,0,0]';
    lam2_opts.degs = lam2_deg; % degree of monomials (not distributed) used to define the operator lam2
    lam2_opts.psatz = 0;
    [prog,lam2] = piesos_sosvar(prog,Zg2,lam2_opts);
    dV_up = -dV - 2*k*Vx - lam2*g;
else
    dV_up = -dV - 2*k*Vx;
end

% Enforce upper bound by defining a SOS distributed monomial.
disp(" --- Enforcing negativity of the derivative ---")

% Declare W3 >= 0
% dV_up has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,dV_up.varname);
ZQ_degmat = unique(floor(dV_up.degmat(:,idcs_x)./2),'rows');
ZQ = dmonomials(x,ZQ_degmat);
Q3_opts.degs = q3_deg; % degree of monomials (not distributed) used to define the operator W3
Q3_opts.exclude = [1,0,0]';
Q3_opts.psatz = 1;      % multipliers (b-s), (s-a)
[prog,W3,Q3mat,ZQ3op] = piesos_sosvar(prog,ZQ,Q3_opts);

% Enforce dV = -W <= 0
disp("  --  enforcing upper bound equality on derivative")
prog = piesos_eq(prog,dV_up-W3);


%% 10. Solve the optimization program
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
    gam_sol = piesos_getsol(prog_sol,gam);
end
