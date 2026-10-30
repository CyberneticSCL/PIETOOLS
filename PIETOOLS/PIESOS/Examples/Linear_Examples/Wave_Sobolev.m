%
% Wave Equation
%
% PDE: x_tt = x_ss;
% BCs: x(t,0) = x(t,1) = 0;
%
% To model, we introduce state components u1=x and u2=x_{t}, to get
%
% PDE: u1_{t} = u2;
%      u2_{t} = u1_{ss}
% BCs: u1(t,0) = u1(t,1) = 0;
%
% We also declare the output z = (u1_s,u2), so that the PIE has the
% output PIE.gz = R*x = z.
%
% This script tests stability of the PIE, by constructing a distributed
% SOS Lyapunov functional V(x) = <Zop*Z(x),Pmat*Zop*Z(Tx)>, with
% Z(x) = [x ... x^d]' and Zop a basis of multiplier and integral
% operators (if use_L2_LF=true, V is declared in terms of the PDE state
% u = Tx instead: V(x) = <Zop*Z(u),Pmat*Zop*Z(u)>), and enforcing
%   V >= eppos*||z||^2;      (lower bound)
%   V <= gam*||z||^2;        (upper bound, gam is minimized)
%   dV <= 0;                 (derivative along the PIE)
% by declaring nonnegative distributed polynomial functionals W1, W2 and
% W3 with V - eppos*||z||^2 = W1, gam*||z||^2 - V = W2 and -dV = W3. With
% use_L2_pos=true or use_L2_bnd=true, ||z|| is replaced by the L2 norm
% ||u|| of the PDE state in the lower or upper bound, respectively; this
% L2 upper bound is infeasible for the wave equation.
%
% The user may change the parameters in "Script parameters": the degrees
% d and pdeg of the functional, the flags use_L2_LF, use_L2_pos and
% use_L2_bnd, the coefficient eppos, and the degrees q1_deg, q2_deg and
% q3_deg of the monomials of W1, W2 and W3.
%
% DJ, 10/15/2026: Initial coding;
% DJ, 10/29/2026: Use the regulated output polynomial PIE.gz, instead of PIE.h;

clear;  clear stateNameGenerator

% Declare the nonlinear PDE
pvar   s t s_dum
dom = [0,1];
u1 = pde_var(s,dom);
u2 = pde_var(s,dom);

z = pde_var('out',2,s,dom);  % output z = (u1_s,u2), to obtain Rop: z = Rx 
PDE = [diff(u1,t)==u2;
       diff(u2,t)==diff(u1,s,2);
       subs(u1,s,0)==0; subs(u1,s,1)==0;
       %subs(u2,s,0)==0; subs(u2,s,1)==0;           % if uncommented, set u2 = pde_var(s,dom,2)
       z==[diff(u1,s);u2]];  

% Convert to a PIE
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent'); 
Top = PIE.T;
f = PIE.f;
x = PIE.x;                                                        

% Script parameters
d = 1;               % degree of LF distributed monomial basis (will be doubled in LF).
pdeg = 2;            % degree of Zs monomials of positive P operator.  
use_L2_LF = false;   % Express the Lyapunov function in terms of an L2 norm of the PDE state
use_L2_pos = false;   % Enforce a lower bound on the LF in terms of the L2 norm of the PDE state
use_L2_bnd = false;  % Upper bound on the LF in terms of the output norm ||z|| (L2 bound infeasible for the wave equation) 
eppos = 0.1; % coefficient in LF lower bound condition.
q1_deg = 2;                                                              
q2_deg = 2; % degree of monomials (not distributed) used to define the operator W2 
q3_deg = 4; % degree of monomials (not distributed) used to define the operator W3 




%% 2. Setting up PIESOS program

% Initialize PIESOS program structure.
dpvar gam
prog = piesos_program(x,gam); % gam

% % Define as a minimization problem on the upper bound (gamma) of the LF 
% % over the domain L_{2,R}.
prog = piesos_setobj(prog,gam);  

%% 3. Construct LF and derivative
% V(v) = <Z_d(v), Pop*Z_d(Tv)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Declare the Lyapunov functional.
disp(" --- Declaring the distributed SOS Lyapunov functional ---")

if use_L2_LF
    Popts.degs = struct('in',pdeg,'out',0);
    Popts.exclude = [0;0;0];
    [prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pde',(1:d)',Popts);
    Tx = Top*x;
else
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
end


%% 4. Define the lower bound on the LF
if use_L2_pos % LF lower bound dependent on ball.
    V_low = Vx - eppos*innerprod(Tx,Tx);
else
    Rx = PIE.gz;
    V_low = Vx - eppos*innerprod(Rx,Rx);
end

% Enforce lower bound by defining a SOS distributed monomial
disp(" --- Enforcing the bound on the Lyapunov functional ---")

% Declare W1 >= 0
% V_low has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,V_low.varname);
Z_bnd_degmat = unique(floor(V_low.degmat(:,idcs_x)./2),'rows'); % This is the degree of W1 
Z_bnd_degmat = Z_bnd_degmat(sum(Z_bnd_degmat,2)>0,:);    % no constant monomials
Z_bnd = dmonomials(x,Z_bnd_degmat); 
Q1_opts.degs = q1_deg; % degree of monomials (not distributed) used to define the operator W1
Q1_opts.exclude = [1,0,0]';
Q1_opts.psatz = 1;      % multipliers (b-s), (s-a) 
[prog,W1,Q1mat,ZQ1op] = piesos_sosvar(prog,Z_bnd,Q1_opts);

disp("  --  enforcing lower bound equality")
prog = piesos_eq(prog,V_low-W1);


%% 5. Define the upper bound on the LF over the specified ball and enforce constraint.
if use_L2_bnd % LF upper bound dependent on ball.
    V_up = gam*innerprod(Tx,Tx)-Vx;
else
    V_up = gam*innerprod(Rx,Rx)-Vx;
end

% Enforce upper bound by defining a SOS distributed monomial
disp(" --- Enforcing the bound on the Lyapunov functional ---")

% Declare W2 >= 0
% V_up has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,V_up.varname);
Z_bnd_degmat = unique(floor(V_up.degmat(:,idcs_x)./2),'rows'); 
Z_bnd_degmat = Z_bnd_degmat(sum(Z_bnd_degmat,2)>0,:);    % no constant monomials
Z_bnd = dmonomials(x,Z_bnd_degmat);                                         
Q2_opts.degs = q2_deg; % degree of monomials (not distributed) used to define the operator W2
Q2_opts.exclude = [1,0,0]';
Q2_opts.psatz = 1;      % multipliers (b-s), (s-a)                          
[prog,W2,Q2mat,ZQ2op] = piesos_sosvar(prog,Z_bnd,Q2_opts); 

disp("  --  enforcing upper bound equality")
prog = piesos_eq(prog,V_up-W2);


%% 6. Define the upper bound on the LF derivative over the specified ball and enforce constraint.

% Enforce upper bound by defining a SOS distributed monomial.
disp(" --- Enforcing negativity of the derivative ---")

% Declare W3 >= 0
% dV has unused variables besides the PIE state x, keep only the columns of x
[~,idcs_x] = ismember(x.varname,dV.varname);    
ZQ_degmat = unique(floor(dV.degmat(:,idcs_x)./2),'rows');
ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);    % no constant monomials 
ZQ = dmonomials(x,ZQ_degmat);  
Q3_opts.degs = q3_deg; % degree of monomials (not distributed) used to define the operator W3 
Q3_opts.exclude = [0,0,0]';
Q3_opts.psatz = 1;      % multipliers (b-s), (s-a)   
[prog,W3] = piesos_sosvar(prog,ZQ,Q3_opts);

% Enforce dV = W3 >= 0
disp("  --  enforcing upper bound equality on derivative")
prog = piesos_eq(prog,-dV-W3);


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
