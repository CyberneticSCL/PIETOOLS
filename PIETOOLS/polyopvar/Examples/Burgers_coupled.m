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
pvar s s_dum t
dom = [0,0.5];
r = 0;

% Declare a PDE and convert to PIE to get the relevant maps from PIE state
% to PDE state
x1 = pde_var(s,dom);       z1 = pde_var('out',s,dom);
x2 = pde_var(s,dom);       z2 = pde_var('out',s,dom);
Dyn = [diff(x1,t)==diff(x1,s,2);
       diff(x2,t)==diff(x2,s,2);
       z1==diff(x1,s);
       z2==diff(x2,s)];
% BCs = [subs(x1,s,dom(1))==0;       subs(x2,s,dom(2))==0;
%        subs(x1,s,dom(2))==0;       subs(x2,s,dom(1))==0];
BCs = [subs(x1,s,dom(1))==0;       subs(x2,s,dom(2))==0;
       subs(x1,s,dom(2))==subs(x2,s,dom(1));
       subs(diff(x1,s),s,dom(2))==subs(diff(x2,s),s,dom(1))];
PIE = convert([Dyn;BCs]);
T = PIE.T;      % [x1;x2] = T*[x1_{ss}; x2_{ss}];
R = 0*PIE.C1;     % [x1_{s};x2_{s}] = T*[x1_{ss}; x2_{ss}];
T11 = dopvar2ndopvar(T(1,1));
T12 = dopvar2ndopvar(T(1,2));
T21 = dopvar2ndopvar(T(2,1));
T22 = dopvar2ndopvar(T(2,2));
R11 = dopvar2ndopvar(R(1,1));
R12 = dopvar2ndopvar(R(1,2));
R21 = dopvar2ndopvar(R(2,1));
R22 = dopvar2ndopvar(R(2,2));

% Declare identity nopvar for identity operators present in polyopvar 
% representation of PIE.
id = nopvar();
id.C = {1, 0, 0};
id.deg = 0;
id.dom = dom;
id.vars = [s, s_dum];

% Declare the right-hand side of the PIE associated with the PDE
%   x1_{t} = x1_{ss} + r*x1 - x1*x1_{s} 
%          = id*x1_f + r*T11*x1_f + r*T12*x2_f - (T11*x1_f)*(R11*x1_f)
%                       - (T12*x2_f)*(R11*x1_f) - (T11*x1_f)*(R12*x2_f)
%                           - (T12*x2_f)*(R12*x2_f)
f = polyopvar();
f.varname = {'x1';'x2'};
f.pvarname = {'s'};
f.dom = dom;
f.varmat = [1;1];
f.degmat = [1,0; 0,1; 2,0; 1,1; 0,2];
C = tensopvar();
C.ops = {id+r*T11, r*T12, {-T11,R11}, {-T11,R12; R11,-T12}, {-T12,R12};
       r*T21, id+r*T22, {-T21,R21}, {-T21,R22; R21,-T22}, {-T22,R22}};
f.C = C;


% Declare the PIE as a struct
PIE = struct();
PIE.T = dopvar2ndopvar(T);
PIE.f = f;
%PIE.vars = [T.var1, T.var2];
%PIE.dom = dom;


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = lpiprogram(s,t,dom);

% Declare monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
vartab = polyopvar({'x1';'x2'},s,dom); % Z_d(v).
Z = dmonomials(vartab,(1:d));

% Declare PSD operator acting on degree d monomial basis and add variables to LPI program.
% opts = ;
pdegs = 0; % maximal monomial degree of Zop. 
Popts.exclude = [0;1;1];
[prog,Pmat,Zop] = soslpivar(prog,Z,pdegs,Popts);
Vx = quad2lin(Pmat,Zop,Z);          

% Add PD constraint to LPI program.
eppos = 1e-4;
Vx.C.ops{1}.params(1,end) = Vx.C.ops{1}.params(end) + eppos;  % can be done more elegantly once polyopvar is better developed
Vx.C.ops{3}.params(1,end) = Vx.C.ops{3}.params(end) + eppos; 

% Evaluate the derivative of V along the PIE
dV = Liediff(Vx,PIE);

% Declare a nonnegative distributed polynomial functional W
qdegs = pdegs+2;
ZQ_degmat = unique(floor(dV.degmat./2),'rows');
ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
ZQ = vartab;
ZQ.degmat = ZQ_degmat;
Q_opts.exclude = [1,0,0]';
[prog,Qmat,ZQop] = soslpivar(prog,ZQ,qdegs,Q_opts);
W = quad2lin(Qmat,ZQop,ZQ);

% Enforce dV = -W <= 0
prog = soslpi_eq(prog,dV+W);

% Solve the optimization program
prog_sol = lpisolve(prog);

