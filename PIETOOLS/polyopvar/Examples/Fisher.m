%
% Burgers-Fishers eq.
%
% PDE: u_t = u_ss + alp*u -u*u_{s} - bet*u^2;
% BCs: u(t,0) = u(t,1) = 0;
%
% should be globally stable when 
%   (1/pi^2)*exp(-3*bet)*(alp+(9/2)*bet^2) < 1
% or perhaps, using a weighted Poincare inequality,
%   (alp+(9/4)*bet^2) <1
% using
%   V(u) = int_{0}^{1} exp(-3*bet*s)*u(s)^2 ds
%
% fundamental state: v = u_ss \in L_2[0,1];
%
%
% PIE: (Tv_t) = f(v) = v + (alp*Tv) -Tv*Rv - bet*(Tv)^2;
%
%

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
clear;  clear stateNameGenerator
pvar s t
dom = [0,1];
alp = 5;
bet = -1;
R = 0.5;          % Radius of ball in which to test stability

% Declare the nonlinear PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2)+alp*x-bet*x^2;
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0];

% Convert to a PIE
PIE = convert(PDE);
Top = PIE.T;
Top_opvar = ndopvar2dopvar(Top);
Rop = dopvar2ndopvar(diff(Top_opvar,Top_opvar.var1,1,'pure'));
f = PIE.f;
x = polyopvar(f.varname,s,dom);


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% % Declare settings
% Set distributed monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
Z = dmonomials(x,(1:d)); % Z_d(v).
% Set distributed monomial basis for psatz multiplier
d_psatz = 1;
Zg = dmonomials(x,(1:d_psatz));
% Monomial degree in independent variables, used to parameter Pop
pdeg = 3;
% Coercivity of Pop
eppos = 1e-2;
% Enforce positivity using piesos_ineq
use_ineq1 = true;
use_ineq2 = true;


% % Initialize PIESOS program structure.
dpvar gam
prog = piesos_program(x,gam);
% Minimize the upper bound gamma on the Lyapunov function
%prog = lpisetobj(prog,gam);
%prog = piesos_ineq(prog,gam);


% % Declare the Lyapunov functional and its derivative,
% % using one of two options
use_opt = 2;
if use_opt==1
    % OPTION 1: Use pie_sosvar
    % This function supports terms of the form <u,u> = int_0^1 u(s)^2 ds
    % but does not support terms of the form <u^2,u^2> = int_0^1 u(s)^4 ds

    % % Declare the Lyapunov functioncal
    disp(" --- Declaring the distributed SOS Lyapunov functional ---")
    % NOTE: This V is expressed in terms of the PDE state
    P_opts.deg = pdeg; % maximal monomial degree of Zop. 
    P_opts.exclude = [0;0;0];
    P_opts.sep = 1;
    [prog,Vx,Pmat,Zop] = piesos_sosvar(prog,Z,P_opts);
    
    % Ensure strict positivity of the Lyapunov functional.
    Vx = Vx + eppos*innerprod(x,x);
                                                    
    % % Compute the derivative of V along the PIE
    disp(" --- Computing the derivative along the PIE ---")
    % NOTE: dV is expressed in terms of the PIE state
    dV = Liediff(Vx,PIE); % output is in polyopvar 

else
    % % OPTION 2: Declare pseudo-qaudratic function directly in terms of
    % ``naive'' monomial basis, [u;u^2;u^3], rather than distributed
    % monomial basis, [u; (u o u); (u o u o u)]:
    %       V(u) = <[u; (u o u)], P [u; u o u]>_{L2}
    %   where P is just a multiplier operator. Then, e.g.
    %       d/dt <(u o u), P (u o u)>_{L2}
    %           = 2*<(T v o Tv), P(Tv o f(v))>_{L2} + 2*<(T v o Tv), P(f(v) o Tv)>_{L2}
    
    % % First, parameterize an SOS multiplier operator
    disp(" --- Declaring the distributed SOS Lyapunov functional ---")
    Zs = cell(d,1);
    Pdim_arr = ones(d,1);
    for i=1:d
        Zs{i} = s.^(0:max(pdeg*i,0));
    end
    [prog,Pcell] = sosquadvar(prog,Zs,Zs,Pdim_arr,Pdim_arr,'pos');
    % Ensure strict positivity of the Lyapunov functional
    Pcell{1} = Pcell{1} + eppos*eye(size(Pcell{1}));
    % Construct the function V in terms of the PIE state
    Tx = Top*x;
    Vx = innerprod(Tx,Tx,Pcell{1});
    if d>=2
        TTx = Tx*Tx;
        Vx = Vx + 2*innerprod(Tx,TTx,Pcell{1,2});
        Vx = Vx + innerprod(TTx,TTx,Pcell{2,2});
    end
    if d>=3
        TTTx = Tx*Tx*Tx;
        Vx = Vx + 2*innerprod(Tx,TTTx,Pcell{1,3});
        Vx = Vx + 2*innerprod(TTx,TTTx,Pcell{2,3});
        Vx = Vx + innerprod(TTTx,TTTx,Pcell{3,3});
    end

    % % Then, compute the derivative along the PIE
    disp(" --- Computing the derivative along the PIE ---")
    dV = 2*innerprod(Tx,f,Pcell{1,1});
    if d>=2
        TTx = Tx*Tx;
        Tfx = Tx*f;     fTx = f*Tx;
        dV12 = 2*innerprod(Tx,Tfx,Pcell{1,2}) + 2*innerprod(Tx,fTx,Pcell{1,2});
        dV22 = 2*innerprod(TTx,Tfx,Pcell{2,2}) + 2*innerprod(TTx,fTx,Pcell{2,2});
        dV = dV + dV12 + dV22;
    end
    if d>=3
        TTTx = Tx*Tx*Tx;
        TTfx = TTx*f;     TfTx = Tx*f*Tx;   fTTx = f*TTx;
        dV13 = 2*(innerprod(Tx,TTfx,Pcell{1,3}) + innerprod(Tx,TfTx,Pcell{1,3}) + innerprod(Tx,fTTx,Pcell{1,3}));
        dV23 = 2*(innerprod(TTx,TTfx,Pcell{2,3}) + innerprod(TTx,TfTx,Pcell{2,3}) + innerprod(TTx,fTTx,Pcell{2,3}));
        dV33 = 2*(innerprod(TTTx,TTfx,Pcell{3,3}) + innerprod(TTTx,TfTx,Pcell{3,3}) + innerprod(TTTx,fTTx,Pcell{3,3}));
        dV = dV + dV13 + dV23 +dV33;
    end
end

% % Add the psatz term, to allow dV>0 when g = R-<Tx,Tx> <=0
if R>0
    disp(" --- Declaring the psatz term ---")
    % Declare the function defining the semialgebraic set
    g = R^2-innerprod(Top*x,Top*x);                             % sphere in L2 norm
    %g = R^2-(innerprod(Top*x,Top*x) +innerprod(Rop*x,Rop*x));   % sphere in Sobolev norm, doesn't really work...
    
    % Declare the multiplier for the upper bound
    lam1_opts.exclude = [1,0,0]';
    lam1_opts.deg = pdeg;
    lam1_opts.psatz = 0;
    [prog,lam1] = piesos_sosvar(prog,Zg,lam1_opts);
    V_bnd = gam*innerprod(Tx,Tx)-Vx-lam1*g;

    % Declare the multiplier for the derivative
    lam2_opts.exclude = [1,0,0]';
    lam2_opts.deg = pdeg;
    lam2_opts.psatz = 0;
    [prog,lam2] = piesos_sosvar(prog,Zg,lam2_opts);
    Wgg = lam2*g;
    dV_g = dV+Wgg;
else
    % Global stability test
    dV_g = dV;
end


% % Enforce Mx-Vx >= 0 using distributed SOS
disp(" --- Enforcing the bound on the Lyapunov functional ---")
if use_ineq1
    disp("  --  using piesos_ineq")
    ineq1_opts.psatz = 0;
    [prog,W1,Q1mat,ZQ1op] = piesos_ineq(prog,V_bnd,ineq1_opts);
else
    % Declare W >= 0
    disp("  --  manually declaring a distributed SOS functional")
    % Determine monomials to balance those of dV
    Z_bnd_degmat = unique(floor(dV_g.degmat./2),'rows');
    Z_bnd = polyopvar(f.varname,s,dom);
    Z_bnd.degmat = Z_bnd_degmat;
    % Declare distributed SOS functional in terms of these monomials
    Q1_opts.deg = pdeg;
    Q1_opts.exclude = [1,0,0]';
    Q1_opts.psatz = 0:1;
    [prog,W1,Q1mat,ZQ1op] = piesos_sosvar(prog,Z_bnd,Q1_opts);

    % Enforce gam*||u||^2-V = W >= 0
    disp("  --  enforcing equality")
    prog = piesos_eq(prog,V_bnd-W1);
end


% % Enforce dV <= 0 using distributed SOS
disp(" --- Enforcing negativity of the derivative ---")
if use_ineq2
    disp("  --  using piesos_ineq")
    ineq2_opts.psatz = 0;
    [prog,W2,Q2mat,ZQ2op] = piesos_ineq(prog,-dV_g,ineq2_opts);
else
    % Declare W >= 0
    disp("  --  manually declaring a distributed SOS functional")
    % Determine monomials to balance those of dV
    ZQ_degmat = unique(floor(dV_g.degmat./2),'rows');
    ZQ = polyopvar(f.varname,s,dom);
    ZQ.degmat = ZQ_degmat;
    % Declare distributed SOS functional in terms of these monomials
    Q2_opts.deg = pdeg+1;
    Q2_opts.exclude = [1,0,0]';
    Q2_opts.psatz = 0:1;
    [prog,W2,Q2mat,ZQ2op] = piesos_sosvar(prog,ZQ,Q2_opts);
    
    % Enforce dV = -W <= 0
    disp("  --  enforcing equality")
    prog = piesos_eq(prog,dV_g+W2);
end


% % Solve the optimization program
disp(" --- Solving the optimization program ---")
sol_opts.simplify = true;
prog_sol = lpisolve(prog,sol_opts);


% % Extract the solution
sol_info = prog_sol.solinfo.info;
if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
    disp("PIESOS program was not solved.")
else
    disp("PIESOS program was successfully solved.")
    Vsol = piesos_getsol(prog_sol,Vx);
    gam_sol = piesos_getsol(prog_sol,gam);
% Pmat_sol = piesos_getsol(prog_sol,Pcell)
% Wg_sol = piesos_getsol(prog_sol,Wg);
% Wnew_sol = piesos_getsol(prog_sol,Wnew);
% dV_sol = piesos_getsol(prog_sol,dV);
end


% % REMARKS
% If we have feasibility on ball of radius R in L2 norm, then this
% certifies stability on any closed level set of V contained in that ball. 
% Consider the gamma-level set of V, {u | V(u) <= gamma}.
% Since V(u)>=eppos*||u||_{L2}, we have ||u||_{L2} <= gamma/eppos for all u
% in the level set. Thus, choosing gamma=eppos*R, the gamma level set of V
% is contained entirely in the ball. Now, we just need a ball that is
% contained in this level set. Suppose that we have
% V(u) <= sqrt{||Pop||_{op}} ||u||_{L2}. Thus, the ball 
% ||u||_{L2} <= gamma/sqrt{||Pop||_op} is contained entirely in the gamma 
% level set of V. Setting gamma = eppos*R, this suggests we can verify
% stability on the ball
%   ||u||_{L2} <= R*eppos/sqrt{||Pop||_{op}}
%
% In our case, however, we have
%   V(u) = <[u;u^2],Pop*[u;u^2]> <= ||Pop||*||[u;u^2]||_{L2}^2
% Thus, we need to bound ||u^2||_{L2}^2 by ||u||_{L2}...