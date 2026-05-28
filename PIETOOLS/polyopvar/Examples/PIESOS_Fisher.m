%
% Fishers eq.
%
% PDE: u_t = u_ss + alp*u - bet*u^2;
% BCs: u(t,0) = u(t,1) = 0;
%
% For alp=5 and bet=-1;
% Locally stable on ball in L2 norm of radius R ~ 4.405
%
%

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
PIE = convert(PDE);
Top = PIE.T;
Top_opvar = ndopvar2dopvar(Top);
Rop = dopvar2ndopvar(diff(Top_opvar,Top_opvar.var1,1,'pure'));
f = PIE.f;
x = f.vartab;


%%%% 2. Test stability using PIESOS

% % Declare settings
% Set distributed monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
Z = dmonomials(x,(1:d)); % Z_d(v).
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
k = 0;      % Rate of exponential decay
% Enforce local decay on semialgebraic set where g>=0
R = 4.01;   % feasible up to R~4.0479 
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
% % In order to declare this functional directly in terms of the
% % fundamental state, rather than the PDE state, we use a manual
% % declaration:
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
    %Zs{i} = 1;
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


% % Add the psatz term, to allow dV>0 when g = R-<Tx,Tx> <=0
if R>0
    disp(" --- Declaring the psatz term ---")
    % Declare the multiplier for the upper bound
    if use_bnd
        if d_psatz1>1
            lam1_opts.exclude = [1,0,0]';
            lam1_opts.deg = pdeg;
            lam1_opts.psatz = 0;
            [prog,lam1] = piesos_sosvar(prog,Zg1,lam1_opts);
            V_bnd = gam*innerprod(Tx,Tx)-Vx-lam1*g;  %(gam*||Tx||^2 - Vx) >= lam1*g >= 0
        else
            V_bnd = gam*innerprod(Tx,Tx)-Vx;
        end
    end

    % Declare the multiplier for the derivative
    lam2_opts.exclude = [1,0,0]';
    lam2_opts.deg = pdeg;
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
        Z_bnd_degmat = unique(floor(V_bnd.degmat./2),'rows');
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
    ZQ_degmat = unique(floor(dV_g.degmat./2),'rows');
    ZQ = polyopvar(f.varname,s,dom);
    ZQ.degmat = ZQ_degmat;
    % Declare distributed SOS functional in terms of these monomials
    qdeg = pdeg+1;
    if d_psatz2==2
        Q2_opts.deg = {[qdeg,qdeg]; [qdeg,round(qdeg/2)]; [qdeg-2,round(qdeg/4)]};
    else
        Q2_opts.deg = qdeg;
    end    
    % Q2_opts.deg = pdeg+1;
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
sol_opts.solver = 'mosek';
prog_sol = lpisolve(prog,sol_opts);


% % Extract the solution
sol_info = prog_sol.solinfo.info;
if sol_info.pinf || sol_info.dinf || sol_info.numerr || abs(sol_info.feasratio-1)>0.1
    disp("PIESOS program was not solved.")
else
    disp("PIESOS program was successfully solved.")
    Vsol = piesos_getsol(prog_sol,Vx);
    gam_sol = piesos_getsol(prog_sol,gam);
end