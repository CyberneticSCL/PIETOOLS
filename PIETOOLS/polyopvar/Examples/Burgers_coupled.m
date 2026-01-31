%
% Coupled Burgers equations
%
% PDE: u1_t = u1_ss + r*u1 - u1*u1_s;
%      u2_t = u2_ss + r*u2 - u2*u2_s;
% BCs: u1(t,0) = 0;     u1(t,0.5) = u2(t,0);
%      u2(t,0.5) = 0;   u1_{s}(t,0.5) = u2_{s}(t,0);
%
% fundamental state: v = [v1;v2] = [u1_ss; u2_ss] \in L_2[0,1];
%
%
% PIE: (T(1,:)*v_t) = f1(v) = v1 + (r*T(1,:)*v) - (T(1,:)*v)(R(1,:)*v);
%      (T(2,:)*v_t) = f2(v) = v2 + (r*T(2,:)*v) - (T(2,:)*v)(R(2,:)*v);
%

%%%% 1. Modelling PIE.

% declare spatial variables, domain, and r.
pvar s s_dum t
dom = [0,0.5];
r = 1;

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
PIE_lin = convert([Dyn;BCs]);
T = PIE_lin.T;      % [x1;x2] = T*[x1_{ss}; x2_{ss}];
R = PIE_lin.C1;     % [x1_{s};x2_{s}] = T*[x1_{ss}; x2_{ss}];
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
% f.degmat = [1,0; 0,1];
% C = tensopvar();
% C.ops = {id+r*T11, r*T12;
%        r*T21, id+r*T22};
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
prog = lpiprogram([s,s_dum],dom);

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
is_diag_term = sum(Vx.degmat,2)==2 & max(Vx.degmat,[],2)==2;
for i=find(is_diag_term')
    Vx.C.ops{i}.params(end) = Vx.C.ops{i}.params(end) + eppos;  % can be done more elegantly once polyopvar is better developed
end

% Evaluate the derivative of V along the PIE
dV = Liediff(Vx,PIE);

% Declare a nonnegative distributed polynomial functional W
qdegs = pdegs+2;
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






% Zop_opvar = blkdiag(ndopvar2dopvar(Zop.ops{1,1}), ndopvar2dopvar(Zop.ops{2,2}));
% Pop = Zop_opvar'*Pmat*Zop_opvar;
% % 
% % tst = subs(Vx.C.ops{1}.params(1),Vx.C.ops{1}.vars,PIE.vars');
% % Check that the integral terms match
% Vparam22 = subs(Vx.C.ops{1}.params(1),Vx.C.ops{1}.vars,[Pop.var1;Pop.var2]);
% Vparam22_alt = 2*Pop.R.R1(2,2);
% diff_22 = Vparam22 - Vparam22_alt;
% Vparam12 = subs(Vx.C.ops{2}.params(1),Vx.C.ops{2}.vars,[Pop.var1;Pop.var2]);
% Vparam12_alt = 2*Pop.R.R1(1,2);
% diff_12 = Vparam12 - Vparam12_alt;
% Vparam21 = subs(Vx.C.ops{2}.params(2),Vx.C.ops{2}.vars,[Pop.var1;Pop.var2]);
% Vparam21_alt = 2*Pop.R.R2(1,2);
% diff_21 = Vparam21 - Vparam21_alt;
% Vparam11 = subs(Vx.C.ops{3}.params(1),Vx.C.ops{3}.vars,[Pop.var1;Pop.var2]);
% Vparam11_alt = 2*Pop.R.R1(1,1);
% diff_11 = Vparam11 - Vparam11_alt;
% 
% % Check that the multiplier terms match
% Vparam22m = subs(Vx.C.ops{1}.params(end),Vx.C.ops{1}.vars(1),T.var1);
% Vparam22m_alt = Pop.R.R0(2,2);
% diff_22m = Vparam22m - Vparam22m_alt;
% Vparam12m = subs(Vx.C.ops{2}.params(end),Vx.C.ops{2}.vars(1),T.var1);
% Vparam12m_alt = 2*Pop.R.R0(1,2);
% diff_12m = Vparam12m - Vparam12m_alt;
% Vparam11m = subs(Vx.C.ops{3}.params(end),Vx.C.ops{3}.vars(1),T.var1);
% Vparam11m_alt = Pop.R.R0(1,1);
% diff_11m = Vparam11m - Vparam11m_alt;
% 
% 
% % 
% eppos = 1e-4;
% is_diag_term = sum(Vx.degmat,2)==2 & max(Vx.degmat,[],2)==2;
% for i=find(is_diag_term')
% Vx.C.ops{i}.params(end) = Vx.C.ops{i}.params(end) + eppos;  % can be done more elegantly once polyopvar is better developed
% end
% Pop.R.R0 = Pop.R.R0 + eppos*eye(size(Pop));
% 
% dV = Liediff(Vx,PIE);
% APT = PIE_lin.A'*Pop*PIE_lin.T + PIE_lin.T'*Pop*PIE_lin.A;
% % Check that the integral terms match
% dVparam22 = subs(dV.C.ops{1}.params(1),dV.C.ops{1}.vars,[APT.var1;APT.var2]);
% dVparam22_alt = 2*APT.R.R1(2,2);
% diff_22 = dVparam22 - dVparam22_alt;
% dVparam12 = subs(dV.C.ops{2}.params(1),dV.C.ops{2}.vars,[APT.var1;APT.var2]);
% dVparam12_alt = 2*APT.R.R2(1,2);
% diff_12 = dVparam12 - dVparam12_alt;
% dVparam21 = subs(dV.C.ops{2}.params(2),dV.C.ops{2}.vars,[APT.var1;APT.var2]);
% dVparam21_alt = 2*APT.R.R1(1,2);
% diff_21 = dVparam21 - dVparam21_alt;
% dVparam11 = subs(dV.C.ops{3}.params(1),dV.C.ops{3}.vars,[APT.var1;APT.var2]);
% dVparam11_alt = 2*APT.R.R1(1,1);
% diff_11 = dVparam11 - dVparam11_alt;
% % tst = subs(dV.C.ops{1}.params(1),dV.C.ops{1}.vars,PIE.vars');
% % tst2 = tst - 2*APT.R.R1;
% % max(max(abs(tst2.C)))
% 
% 
% % ZQop_opvar = ndopvar2dopvar(ZQop.ops{1});
% % Qop = ZQop_opvar'*Qmat*ZQop_opvar;
% % tst = subs(W.C.ops{1}.params(1),W.C.ops{1}.vars,PIE.vars');
% % tst2 = tst - 2*Qop.R.R1;
% % max(max(abs(tst2.C)))