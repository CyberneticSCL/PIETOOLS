% polynomial heat eq.
% u_t = u_ss + f(s,u) + [-k*u + k*b(s)*(u-u(s^*))];

%% 1. Modelling PIE.

clear;  clear stateNameGenerator

p=2; % degree of polynomial term in PDE.
d=1; % degree of SOS LF.
NORM = True; % type of norm: True=L2, False=Sobolev.
R = 0.5; % Radius of ball in which to test stability.
L = 1; % Length of PDE domain.

pvar s t
dom = [0,L];

% Declare the nonlinear PDE
x = pde_var(s,dom);
PDE = [diff(x,t)==diff(x,s,2) + x^p;
       subs(x,s,dom(1))==0;  subs(x,s,dom(2))==0]; % Dirchlet BCs.

% Convert to a PIE
PIE = convert(PDE);
Top = PIE.T;
Top_opvar = ndopvar2dopvar(Top);
Rop = dopvar2ndopvar(diff(Top_opvar,Top_opvar.var1,1,'pure'));
f = PIE.f;
x = polyopvar(f.varname,s,dom);


%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
% V(v) = <Z_d(v), Pop*Z_d(T*v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set distributed monomial basis of SOS Lyapunov functional.
Z = dmonomials(x,(1:d)); % Z_d(v).
ZT = dmonomials(Top*x,(1:d)); % Z_d(Tv).

% Set distributed monomial basis for psatz multiplier
d_psatz = 1;
Zg = dmonomials(x,(1:d_psatz));

% Monomial degree in independent variables, used to parameter Pop
pdeg = 3;

% Coercivity of Pop
eppos = 1e-2;

% % Initialize PIESOS program structure.
dpvar gam
prog = piesos_program(x,gam);
% Minimize the upper bound gamma on the Lyapunov function
%prog = lpisetobj(prog,gam);
%prog = piesos_ineq(prog,gam);

 % Declare the Lyapunov functional and its derivative