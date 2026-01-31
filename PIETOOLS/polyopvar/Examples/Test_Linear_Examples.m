%
% Test the example from the PDE library
%%%% 1. Declare the PIE

example_num = 7;
PDE = examples_PDE_library_PIETOOLS_nouinput(example_num);
PIE_lin = convert(PDE);
PIE = PIE2polyPIE(PIE_lin);

xname = PIE.f.varname;
var1 = polynomial(PIE.f.pvarname);
dom = PIE.f.dom;



%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
prog = lpiprogram(PIE.vars,PIE.dom);

% Declare monomial basis of SOS Lyapunov functional.
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
x = polyopvar(xname,var1,dom); % Z_d(v).
Zvec = dmonomials(x,(1:d));

% Declare PSD operator acting on degree d monomial basis and add variables to LPI program.
% opts = ;
pdegs = 4; % maximal monomial degree of Zop. 
Popts.exclude = [0;0;0];
[prog,Pmat,Zop] = soslpivar(prog,Zvec,pdegs,Popts);
Vx = quad2lin(Pmat,Zop,Zvec); % output is in polyopvar

% Add PD constraint to LPI program.
eppos = 1e-4;
is_diag_term = sum(Vx.degmat,2)==2 & max(Vx.degmat,[],2)==2;
for i=find(is_diag_term')
    Vx.C.ops{i}.params(end) = Vx.C.ops{i}.params(end) + eppos;  % can be done more elegantly once polyopvar is better developed
end                                                         % Vx.C is a tensopvar

% Evaluate the derivative of V along the PIE
dV = Liediff(Vx,PIE); % output is in polyopvar 

% Declare a nonnegative distributed polynomial functional W
qdegs = pdegs+2;
ZQ_degmat = unique(floor(dV.degmat./2),'rows');
ZQ_degmat = ZQ_degmat(sum(ZQ_degmat,2)>0,:);
ZQ = polyopvar(xname,var1,dom);
ZQ.degmat = ZQ_degmat;
Q_opts.exclude = [1,0,0]';
[prog,Qmat,ZQop] = soslpivar(prog,ZQ,qdegs,Q_opts);
W = quad2lin(Qmat,ZQop,ZQ);

% Enforce dV = -W <= 0
prog = soslpi_eq(prog,dV+W);

% Solve the optimization program
prog_sol = lpisolve(prog);

