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
r = 0.1;           % should be [−0.6500, 0.7202]

% Declare a PDE and convert to PIE to get the relevant maps from PIE state
% to PDE state
x = pde_var(s,dom);       z1 = pde_var('out',s,dom);
x2 = pde_var(s,dom);       z2 = pde_var('out',s,dom);
Dyn = [diff(x,t)==diff(x,s,4);
       z1==diff(x,s);
       z2==diff(x,s,2)];
% BCs = [subs(x1,s,dom(1))==0;       subs(x2,s,dom(2))==0;
%        subs(x1,s,dom(2))==0;       subs(x2,s,dom(1))==0];
BCs = [subs(x,s,dom(1))==0;       subs(x,s,dom(2))==0;
       subs(diff(x,s),s,dom(1))==0;
       subs(diff(x,s),s,dom(2))==0];
PIE_lin = convert([Dyn;BCs]);
Top = PIE_lin.T;      
Rop = PIE_lin.C1(1,:);     
R2op = PIE_lin.C1(2,:);
T = dopvar2ndopvar(Top);
R = dopvar2ndopvar(Rop);
R2 = dopvar2ndopvar(R2op);


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
f.varname = {'x'};
f.pvarname = {'s'};
f.dom = dom;
f.varmat = 1;
f.degmat = [1;2];
C = tensopvar();
C.ops = {-id-R2, {-T,R; -r*T,T}};
% f.degmat = [1,0; 0,1];
% C = tensopvar();
% C.ops = {id+r*T11, r*T12;
%        r*T21, id+r*T22};
f.C = C;


% Declare the PIE as a struct
PIE = struct();
PIE.T = T;
PIE.f = f;
%PIE.vars = [T.var1, T.var2];
%PIE.dom = dom;


%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = lpiprogram([s,s_dum],dom);

% Declare monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
vartab = polyopvar({'x'},s,dom); % Z_d(v).
Z = dmonomials(vartab,(1:d));

% Declare PSD operator acting on degree d monomial basis and add variables to LPI program.
% opts = ;
pdegs = 4; % maximal monomial degree of Zop. 
Popts.exclude = [0;0;0];
[prog,Pmat,Zop] = soslpivar(prog,Z,pdegs,Popts);
Vx = quad2lin(Pmat,Zop,Z);          

% Add PD constraint to LPI program.
eppos = 1e-4;
is_diag_term = sum(Vx.degmat,2)==2 & max(Vx.degmat,[],2)==2;
for i=find(is_diag_term')
    Vx.C.ops{i}.params(end) = Vx.C.ops{i}.params(end) + eppos;  % can be done more elegantly once polyopvar is better developed
end

% Evaluate the derivative of V along the PIE
dV = Liediff(Vx,PIE);

% Declare a nonnegative distributed polynomial functional W
qdegs = pdegs+4;
ZQ_degmat = unique(floor(dV.degmat./2),'rows');
ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
ZQ = vartab;
ZQ.degmat = ZQ_degmat;
Q_opts.exclude = [1,0,0]';
Q_opts.psatz = 0;
[prog,Qmat,ZQop] = soslpivar(prog,ZQ,qdegs,Q_opts);
W = quad2lin(Qmat,ZQop,ZQ);

% Enforce dV = -W <= 0
prog = soslpi_eq(prog,dV+W);

% Solve the optimization program
prog_sol = lpisolve(prog);


