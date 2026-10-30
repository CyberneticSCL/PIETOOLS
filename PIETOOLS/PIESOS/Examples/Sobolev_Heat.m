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
% This script tests local stability of the PIE, using the energy
% functional V(u) = 0.5*||u||^2 (use_L2_LF=true). It enforces
%   V <= gam*||u||^2;          (upper bound, gam is minimized)
%   dV <= 0;                   (derivative along the PIE)
% on the region G = {u | g(u) >= 0}, using SOS multipliers lam1 and lam2
% of degrees d_psatz1 and d_psatz2 (the conditions hold globally if these
% are 0). The region is the Sobolev ball g = q_sobolev^2 - ||u_s||^2 if
% use_L2_ball=false, and the L2 ball g = q_sobolev^2 - ||u||^2 otherwise.
% If use_L2_bnd=false, the norm ||u|| in the upper bound is replaced by
% ||u_s||.
%
% The user may change the coefficients L, a and k of the PDE and the
% parameters in "Script parameters": the radius q of the region, the flags
% use_L2_ball and use_L2_bnd, the degrees d_psatz1 and d_psatz2, and the
% degrees q2_deg, q3_deg, lam1_deg and lam2_deg of the monomials of the
% SOS variables W2, W3 and multipliers lam1, lam2.
%
% Remark: only use_L2_LF=true is supported. For use_L2_LF=false, pdeg is
% not set and the derivative dV is not computed.
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


% Script parameters
q = 1.65*a*pi^2/L^2;     % For k=2*a*pi^2/L^2, stability can be verified with energy functional when |u|<q=a*pi^2/L^2
q_sobolev = q/L;     % by Sobolev inequality, |u|^2 <= L*||u_{x}||_{L2}^2
d = 1;               % degree of LF distributed monomial basis (will be doubled in LF).
% pdeg = 4;            % degree of Zs monomials of positive P operator.

% CARL: CONDITIONS TO REPLICATE EMILIA'S RESULT
use_L2_LF = true;   % Express the Lyapunov function in terms of an L2 norm of the PDE state
use_L2_ball = false; % Enforce bounds on the LF and derivative locally on an L2 ball (if TRUE) or Sobolev ball (if false) of radius R.

% use_L2_pos = true;   % Enforce a lower bound on the LF in terms of the L2 norm of the PDE state
use_L2_bnd = true;  % Enforce an upper bound on the LF in terms of the L2 norm of the PDE state
d_psatz1 = 1;
d_psatz2 = 1; % degree of distributed monomial in upper bound condition of LF derivative.

% eppos = 1e-2; % coefficient in LF lower bound condition.
% q1_deg = 4;
q2_deg = 4; % degree of monomials (not distributed) used to define the operator W2
q3_deg = 4; % degree of monomials (not distributed) used to define the operator W3
lam1_deg = 4; % degree of monomials (not distributed) used to define the operator lam1
lam2_deg = 4; % degree of monomials (not distributed) used to define the operator lam2

% Stability can be verified with energy functional if Phi<=0, where
% Phi = [q_sobolev-k, k/2;k/2, -a*pi^2/(L^2)];

%% 1. Modelling PIE.

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');
Top = PIE.T;
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

if use_L2_LF
    Tx = Top*x;
    Vx = 0.5*innerprod(Tx,Tx);
else
    % Declare the functional, and its derivative along the PIE, in terms of the
    % PIE state, with
    %   V(x) = <Zop*Z(x), Pmat*Zop*Z(Tx)>,
    % for Pmat constant (not enforcing positivity here), and Zop a basis of
    % multiplier and integral operators with monomials of degree up to pdeg.
    % The symmetry condition V(xL,T*xR) = V(xR,T*xL) is enforced by piesos_LF.
%{
    % Construct monomial basis for Pop=P*Z
    Zmon1 = monomials(s,0:pdeg);
    Zmon2 = monomials([s,s_dum],0:pdeg);
    Zop = opvar();
    Zop.R.R0 = [Zmon1;0*Zmon2;0*Zmon2];
    Zop.R.R1 = [0*Zmon1;Zmon2;0*Zmon2];
    Zop.R.R2 = [0*Zmon1;0*Zmon2;Zmon2];
    Zop.var1 = s;
    Zop.var2 = s_dum;
    Zop.I = dom;
    Z = dopvar2ndopvar(Zop);
    
    % P matrix (not enforcing positivity here)
    [prog,Pcell] = sosquadvar(prog,{1},{1},size(Zop,1),1); % 1's restrict Pcell to matrix (not function)
    
    % Construct LF (up to degree 3).
    Tx = Top*x;
    Zx = Z*x;
    Vx = innerprod(Tx,Zx,Pcell{1}');
%}
    Popts.degs = struct('in',pdeg,'out',0);
    Popts.exclude = [0;0;0];
    [prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pie',(1:d)',Popts);
    Tx = Top*x;
end


%% Enforce the symmetry condition
% The symmetry condition V(xL,T*xR) = V(xR,T*xL) is enforced by piesos_LF.


%% 4. Construct the derivative of the LF along the PIE
disp(" --- Computing the derivative along the PIE ---")
if use_L2_LF
    dV = innerprod(Tx,f);
else
end


%% 5. Define the specified local ball as a semialgebraic set.
disp(" --- Declaring the local ball ---")
if use_L2_ball
    g = q_sobolev^2 - innerprod(Tx,Tx); % L2 ball of radius R.
else
    Rx = PIE.gz;
    g = q_sobolev^2 - innerprod(Rx,Rx); % Sobolev ball of radius R. - innerprod(Tx,Tx) 
end


%% 6. Define the lower bound on the LF (holds globally) and enforce constraint.
% 
% if use_L2_pos % LF lower bound dependent on ball.
%     V_low = Vx - eppos*innerprod(Tx,Tx);
% else
%     Rx = Rop*x;
%     V_low = Vx - eppos*innerprod(Rx,Rx);
% end
% 
% % Enforce lower bound by defining a SOS distributed monomial
% disp(" --- Enforcing the bound on the Lyapunov functional ---")
% 
% % Declare W1 >= 0
% Z_bnd_degmat = unique(floor(V_low.degmat./2),'rows'); % This is the degree of W1
% Z_bnd = polyopvar(f.varname,s,dom);
% Z_bnd.degmat = Z_bnd_degmat;
% Q1_opts.deg = q1_deg; % degree of monomials (not distributed) used to define the operator W1
% Q1_opts.exclude = [1,0,0]';
% Q1_opts.psatz = 0:1;
% [prog,W1,Q1mat,ZQ1op] = piesos_sosvar(prog,Z_bnd,Q1_opts);
% 
% disp("  --  enforcing lower bound equality")
% prog = piesos_eq(prog,V_low-W1);


%% 7. Define the upper bound on the LF over the specified ball and enforce constraint.

% Declare the SOS multiplier, lam1, then define bound.
Zg1 = dmonomials(x,(1:d_psatz1));
if d_psatz1>=1
    lam1_opts.exclude = [1,0,0]';
    lam1_opts.degs = lam1_deg; % degree of monomials (not distributed) used to define the operator lam1
    lam1_opts.psatz = 0;
    [prog,lam1] = piesos_sosvar(prog,Zg1,lam1_opts);

    if use_L2_bnd % LF upper bound dependent on ball.
        V_up = gam*innerprod(Tx,Tx)-Vx-lam1*g;
    else
        V_up = gam*innerprod(Rx,Rx)-Vx-lam1*g;
    end
else
    if use_L2_bnd % LF upper bound dependent on ball.
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
[prog,W2,Q2mat,ZQ2op] = piesos_sosvar(prog,Z_bnd,Q2_opts); % DO WE NEED Q2MAT AND ZQ2OP??????????????????

disp("  --  enforcing upper bound equality")
prog = piesos_eq(prog,V_up-W2);

%% 9. Define the upper bound on the LF derivative over the specified ball and enforce constraint.

% Declare the SOS multiplier, lam2, then define bound.
Zg2 = dmonomials(x,(1:d_psatz2));
if d_psatz2>=1
    lam2_opts.exclude = [1,0,0]';
    lam2_opts.degs = lam2_deg; % degree of monomials (not distributed) used to define the operator lam2
    lam2_opts.psatz = 0;
    [prog,lam2] = piesos_sosvar(prog,Zg2,lam2_opts);
    dV_up = -dV - lam2*g;
else
    dV_up = -dV;
end

% Enforce upper bound by defining a SOS distributed monomial.
disp(" --- Enforcing negativity of the derivative ---")

% Declare W3 >= 0
% dV_up has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,dV_up.varname);
ZQ_degmat = unique(floor(dV_up.degmat(:,idcs_x)./2),'rows'); % CARL: THIS CHOICE OF ZQ_degmat DOES NOT ALLOW
% W3.degmat TO MATCH dV_up.degmat. In this case W3.degmat = [0 1 2] whilst dV_up.degmat = [0 2 3] 
% --> IMPOSSIBLE TO SATISFY EQUALITY CONSTRAINT.
% THIS PROBLEM WAS ONLY TRUE WHEN dV_up = -dV - g;
% ZQ_degmat = unique(ceil(dV_up.degmat./2),'rows'); % --> W3.degmat = [0 1 2 3 4]

ZQ = dmonomials(x,ZQ_degmat);
Q3_opts.degs = q3_deg; % degree of monomials (not distributed) used to define the operator W3
Q3_opts.exclude = [0,0,0]';
Q3_opts.psatz = 1;      % multipliers (b-s), (s-a)
[prog,W3] = piesos_sosvar(prog,ZQ,Q3_opts);

% Enforce dV = W3 >= 0
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
    % gam_sol = piesos_getsol(prog_sol,gam);
end


% Compute the level set
%gam_sol = GetLevelSet(Vsol,g);