%
% Test the stability of an example from the PDE library
%
% This script takes the PDE with number example_num from the library
% `examples_PDE_library_PIETOOLS_nouinput`, converts it to a PIE with
% `convert_PIETOOLS_PDE_nonlinear`, and tests stability of the PIE. It
% declares a distributed SOS Lyapunov functional V(x) =
% <Zop*Z(u),Pmat*Zop*Z(u)> + eppos*||u||^2, with u = T*x the PDE state,
% Z(u) = [u ... u^d]', Zop a basis of multiplier and integral operators,
% and Pmat >= 0, and computes its derivative dV along the PIE with
% `piesos_LF`. It then enforces dV <= 0 using `piesos_ineq`. The option
% `psatz` of `piesos_ineq` is the maximal order q of the multipliers
% (b-s)^(q-j)*(s-a)^j, j = 0,...,q, used to enforce positivity of the
% operators on the interval [a,b] of the spatial domain.
%
% The user may change the example number example_num, the degree d, the
% maximal monomial degree pdegs and options Popts (exclude, pos) of the
% functional, the coefficient eppos, and the option ineq_opts.psatz.
%
% DJ, 10/15/2026: Initial coding;

%%%% 1. Declare the PIE

clear
example_num = 15;
PDE = examples_PDE_library_PIETOOLS_nouinput(example_num);
PIE = convert_PIETOOLS_PDE_nonlinear(PDE,'silent');




%%%% 2. Setting up LPI for stability analysis using a SOS Lyapunov functional
%%%%  V(v) = <Z_d(v), Pop*Z_d(v)>; Z_d(v) = [v ... v^d]'; Pop = Zop^* Pmat Zop.

% Set up LPI program structure.
x = PIE.x;
prog = piesos_program(x);

% Declare a SOS Lyapunov functional V, expressed in terms of the PDE state
% u = Top*x, and compute its derivative dV along the PIE, in terms of x:
%   V(x) = <Zop*Z(u),Pmat*Zop*Z(u)>, Pmat>=0, Z(u) = [u ... u^d]';
d = 1;                  % degree of distributed monomial basis, will be doubled in LF
pdegs = 2; % maximal monomial degree of Zop.
Popts.degs = pdegs;
Popts.exclude = [0;0;0];
Popts.pos = true;
[prog,Vx,dV,Pmat,Zop] = piesos_LF(prog,PIE,'pde',(1:d)',Popts);

% Add PD constraint to LPI program, adding eppos*<u,u> to V, with
% derivative 2*eppos*<u,f(x)>
eppos = 1e-4;
u = PIE.T*x;
Vx = Vx + eppos*innerprod(u,u);
dV = dV + 2*eppos*innerprod(u,PIE.f);

% Enforce dV = -W <= 0, for W a nonnegative distributed polynomial functional
ineq_opts.psatz = 1;    % multipliers (b-s), (s-a), see piesos_quadopvar
[prog,W,Qmat,ZQop] = piesos_ineq(prog,-dV,ineq_opts);

% Solve the optimization program
prog_sol = piesos_solve(prog);





% Zop_opvar = ndopvar2dopvar(Zop.ops{1});
% Pop = Zop_opvar'*Pmat*Zop_opvar;
% 
% tst = subs(Vx.C.ops{1}.params(1),Vx.C.ops{1}.vars,PIE.vars');
% tst2 = tst - 2*Pop.R.R1;
% max(max(abs(tst2.C)))
% tmp = Vx.C.ops{1}.params(2);
% tmp2 = subs(tmp,polynomial({'s_1'}),PIE.vars(1)) - Pop.R.R0;
% max(max(abs(tmp2.C)))
% 
% 
% eppos = 1e-4;
% is_diag_term = sum(Vx.degmat,2)==2 & max(Vx.degmat,[],2)==2;
% for i=find(is_diag_term')
% Vx.C.ops{i}.params(end) = Vx.C.ops{i}.params(end) + eppos;  % can be done more elegantly once polyopvar is better developed
% end
% Pop.R.R0 = Pop.R.R0 + eppos;

% dV = Liediff(Vx,PIE);
% APT = PIE_lin.A'*Pop*PIE_lin.T + PIE_lin.T'*Pop*PIE_lin.A;
% tst = subs(dV.C.ops{1}.params(1),dV.C.ops{1}.vars,PIE.vars');
% tst2 = tst - 2*APT.R.R1;
% max(max(abs(tst2.C)))


% ZQop_opvar = ndopvar2dopvar(ZQop.ops{1});
% Qop = ZQop_opvar'*Qmat*ZQop_opvar;
% tst = subs(W.C.ops{1}.params(1),W.C.ops{1}.vars,PIE.vars');
% tst2 = tst - 2*Qop.R.R1;
% max(max(abs(tst2.C)))
