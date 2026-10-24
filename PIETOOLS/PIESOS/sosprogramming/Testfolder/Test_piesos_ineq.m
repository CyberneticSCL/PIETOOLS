function Test_piesos_ineq()
% TEST_PIESOS_INEQ tests, using assertions, piesos_ineq,
%   [prog,g,Qmat,Zop,Z] = piesos_ineq(prog,f,opts),
% which adds the constraint f>=0 to a PIESOS program for a scalar
% 'polyopvar' functional f, by declaring an SOS functional g
% (piesos_sosvar) and imposing f=g (piesos_eq).
%
% Oracle, independent of the function under test:
%  - feasibility: functionals nonnegative for every state (int x^2,
%    int x^4, int x^2 + 2, int s*x^2 on [0,1], (int x)^2, ...) must give a
%    feasible program, functionals negative for some state (-int x^2,
%    int x^2 - 2, int x^2 - int x^4, int (s-1/2)*x^2, an odd power) an
%    infeasible one; g>=0 by construction, so feasibility certifies f>=0;
%  - values: for f=int x^2, the solved Qmat is PSD and its quadratic form
%    integrated against random p by quadrature gives int p^2;
%  - options: rows of Zop with and without a delta term in f;
%    user-specified 'exclude' and 'degs' must not be overwritten;
%  - uncombined f (case_uncombined_f): duplicate monomial rows with
%    cancelling kernels, and a monomial left with a zero operator, must
%    give the same basis size as f written without redundancy;
%  - basis monomials mixing states of different dependence
%    (case_mixed_basis_warning): warning issued only for them; the
%    conservatism (int int (x1*x2)^2 infeasible if x1, x2 depend on
%    different variables, feasible otherwise) is pinned;
%  - kernel degrees (case_kernel_degrees, case_option_types): f SOS by
%    construction with known degrees and term types of Zop.
%  - state variables that do not appear in any monomial, as returned by
%    piesos_LF (case_unused_states): for the transport equation the
%    program is feasible;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_ineq
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/21/2026: Add case_vector_states;
% DJ (CC), 10/23/2026: case_vector_states: basis with a block per state
%   (degrees per state, no longer per scalar row), check its block sizes;
% DJ (CC), 10/23/2026: Merge Test_piesos_ineq_unused_states into this test,
%   as cases;

% NOTES
% - case_mixed_basis_warning pins a known limitation: a basis monomial
%   mixing states of different dependence is discarded, so f=0 is imposed;
%   the test must be updated when mixed monomials are supported;
% - A cross term of Gram matrix [1 -.75;-.75 1] is used in
%   case_two_states: a perfect square (c=-2) is on the boundary of the
%   cone, where SeDuMi reports both pinf and dinf;
% - case_vector_states: states of size larger than 1 (the states, not their
%   rows, are the variables of the monomials);
% - Programs without a PSD variable fail in SeDuMi: feasible() adds a
%   trivial one;
% - The degree rule for Zop is exact for numeric f; f with decision
%   variables (case_symbolic_f) uses the heuristic d=p=ceil((e(1)-1)/3).

tests ={@case_input_errors, @case_sweep_1var, @case_certificate_values, ...
         @case_nonlocal_kernel, @case_user_options, @case_zero_basis, ...
         @case_constants, @case_psatz, @case_two_states, @case_mixed_dependence, ...
         @case_two_variables, @case_uncombined_f, ...
         @case_mixed_basis_warning, @case_kernel_degrees, @case_option_types, ...
         @case_psatz_degree, @case_symbolic_f, @case_vector_states, ...
         @case_unused_states};
nfail = 0;
for i=1:numel(tests)
    try
        tests{i}();
        fprintf('PASS  %s\n',func2str(tests{i}));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_ineq: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_ineq: all %d cases passed.\n',numel(tests));
end

%% ---- cases ------------------------------------------------------------
function case_unused_states()
% case_unused_states -- tests piesos_ineq for a functional
% with state variables that do not appear in any monomial, as returned by
% piesos_LF (names of the PDE state and of the copies of the PIE state,
% besides the PIE state).
%
% For the transport equation x_t = -x_s, x(0)=0, the Lie derivative of the
% functional V = <u,P*u> + eps*<u,u> along the PIE is nonpositive for P=0
% (dV = -eps*u(1)^2), so the program "find P>=0 such that -dV is a sum of
% squares" is feasible. Infeasibility (pinf=1) is therefore a failure.
%

% NOTES
% - A psd variable (sosposmatr) is added so that SeDuMi's preprocessing
%   does not fail on small programs;
% - Feasible means pinf=dinf=0 and feasratio>=0.5.

    clear stateNameGenerator
    pvar s t
    dom = [0,1];
    u = pde_var(s,dom);
    PDE = [diff(u,t)==-diff(u,s);   subs(u,s,dom(1))==0];
    evalc('PIE = convert_PIETOOLS_PDE_nonlinear(PDE,''silent'');');
    x = PIE.x;
    nfeas = 0;

    opts_list = {struct('degs',2,'exclude',[0;0;0]), struct('degs',0,'exclude',[0;1;1]), struct('degs',struct('in',0,'out',2),'exclude',[0;1;1])};
    for io = 1:numel(opts_list)
        prog = piesos_program(x);
        Popts = opts_list{io};      Popts.pos = true;
        evalc('[prog,Vx,dV] = piesos_LF(prog,PIE,''pde'',1,Popts);');
        % Strict positivity: eps*<u,u>, with derivative 2*eps*<u,f(x)>
        Tx = PIE.T*x;
        Vx = Vx + 1e-2*innerprod(Tx,Tx);            %#ok<NASGU>
        dV = dV + 2e-2*innerprod(Tx,PIE.f);
        evalc('[prog,W] = piesos_ineq(prog,-dV);');
        prog = sosposmatr(prog,20);
        evalc('prog = piesos_solve(prog);');
        info = prog.solinfo.info;
        if info.pinf || info.dinf || info.feasratio<0.5
            error("case_unused_states failed (options %d): program infeasible (pinf=%d, feasratio=%g).",io,info.pinf,info.feasratio)
        end
        nfeas = nfeas+1;
    end

    % The variables of dV from piesos_LF are the PIE state only
    prog = piesos_program(x);
    evalc('[~,~,dV] = piesos_LF(prog,PIE,''pde'',1,struct(''degs'',0,''exclude'',[0;1;1],''pos'',true));');
    if ~isequal(reshape(dV.varname,[],1),reshape(x.varname,[],1))
        error("case_unused_states failed: dV has unused state variables.")
    end

    fprintf('case_unused_states passed (%d programs feasible).\n',nfeas)
end

function [prog,x] = new_prog(varnames,name)
% Program in the distributed state x of the given spatial variables, on [0,1]
if nargin<2
    name = 'x1';
end
x = polyopvar(name,polynomial(varnames(:)),[0,1]);
prog = piesos_program(x);
end

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s".',err.message);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function [is_ok,prog] = feasible(prog)
[prog,~] = sosposmatr(prog,1);
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
is_ok = (info.pinf==0 && info.dinf==0);
end

function check_ineq(f,x,expected,label,opts)
% Declare f>=0 for the state x, and check feasibility against expected
prog = piesos_program(x);
if nargin<5
    [prog,g] = piesos_ineq(prog,f);
else
    [prog,g] = piesos_ineq(prog,f,opts);
end
assert(isa(g,'polyopvar'),'%s: g must be a polyopvar.',label);
is_ok = feasible(prog);
assert(is_ok==expected,'%s: program %s, expected %s.',label, ...
    tern(is_ok,'feasible','infeasible'),tern(expected,'feasible','infeasible'));
end

function s = tern(c,a,b)
if c, s = a; else, s = b; end
end

function V = make_V(x,K,d)
% d-th power of the state, with coefficient functional K, as in Test_piesos_eq
V = x;
for i=2:d
    V = V*x;
end
V.C = tensopmat({K});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
end

function case_input_errors()
pvar s
[prog,x] = new_prog({'s'});
f = innerprod(x,x);
assert_err(@() piesos_ineq(5,f),'struct');
assert_err(@() piesos_ineq(prog,'a'),'polyopvar');
assert_err(@() piesos_ineq(prog,f,5),'struct');
assert_err(@() piesos_ineq(prog,f,struct(),1),'Too many input arguments');
% The state itself is a function, not a functional (output variable s)
assert_err(@() piesos_ineq(prog,x*x),'scalar-valued functional');
assert_err(@() piesos_ineq(prog,[f;f]),'scalar-valued functional');
end

function case_sweep_1var()
pvar s
[~,x] = new_prog({'s'});
f2 = innerprod(x,x);                    % int x^2
f4 = innerprod(x.*x,x.*x);              % int x^4
check_ineq(f2,x,true,'int x^2');
check_ineq(-f2,x,false,'-int x^2');
check_ineq(f4,x,true,'int x^4');
check_ineq(f2+f4,x,true,'int x^2 + int x^4');
check_ineq(3*f2+f4,x,true,'3 int x^2 + int x^4');
% Indefinite: negative for large x, positive for small x
check_ineq(f2-f4,x,false,'int x^2 - int x^4');
check_ineq(-f4,x,false,'-int x^4');
end

function case_certificate_values()
% f = int_0^1 x^2: the solved Qmat is PSD, and g(p) = int z(s)'*Qmat*z(s) ds, with
% Zop*p = z = [p; int_0^s p; int_s^1 p], gives int p^2 (quadrature). f has a delta term
% only, so by default Zop = x and Qmat 1 x 1: all three terms need exclude=false
rng(7);
pvar s
[prog,x] = new_prog({'s'});
[~,~,Q1,Zop1] = piesos_ineq(prog,innerprod(x,x));
assert(isequal(size(Q1{1,1}),[1,1]) && size(Zop1,1)==1,'By default, Zop must be x (1 row), Qmat 1 x 1.');
[prog,g,Qmat,Zop,Z] = piesos_ineq(prog,innerprod(x,x),struct('exclude',false(3,1)));
assert(isa(Qmat,'cell') && isequal(size(Qmat),[1,1]) && isequal(size(Qmat{1,1}),[3,3]),'Qmat must be a 1 x 1 cell (1 monomial) with a 3 x 3 block.');
assert(isa(Zop,'tensopmat') && size(Zop,1)==3,'Zop must be a tensopmat with 3 rows.');
assert(isa(Z,'polyopvar') && isequal(full(Z.Z),1),'Z must be the basis [x].');
[ok,prog] = feasible(prog);
assert(ok,'Program must be feasible.');
P = double(piesos_getsol(prog,Qmat{1,1}));
assert(min(eig((P+P')/2))>-1e-6,'Qmat must be positive semidefinite.');
for t=1:3
    c = randn(1,3);                                 % p(s) = c1+c2*s+c3*s^2
    p = @(t) c(1)+c(2)*t+c(3)*t.^2;
    Ip = @(t) c(1)*t+c(2)*t.^2/2+c(3)*t.^3/3;       % int_0^t p
    z = @(t) [p(t); Ip(t); Ip(1)-Ip(t)];
    val = integral(@(t) arrayfun(@(ti) z(ti)'*P*z(ti),t),0,1,'AbsTol',1e-12,'RelTol',1e-10);
    ref = integral(@(t) p(t).^2,0,1,'AbsTol',1e-12,'RelTol',1e-10);
    assert(abs(val-ref)<1e-5*max(1,abs(ref)),'g(p) must equal int p^2: %g vs %g.',val,ref);
end
end

function case_nonlocal_kernel()
% f = +-(int x)^2 = +-int int x(s)x(t): kernel without delta term, so the
% multiplier term is excluded from Zop (2 rows rather than 3)
pvar s
[prog,x] = new_prog({'s'});
K = intop({polynomial(2)},{[0 3 0 1]},{'s'},[0,1]);
f = make_V(x,K,2);
[~,~,Qmat,Zop] = piesos_ineq(prog,f);
assert(size(Zop,1)==2 && isequal(size(Qmat{1,1}),[2,2]),'Without delta term, Zop must have 2 rows, got %d.',size(Zop,1));
check_ineq(f,x,true,'(int x)^2');
f = make_V(x,intop({polynomial(-2)},{[0 3 0 1]},{'s'},[0,1]),2);
check_ineq(f,x,false,'-(int x)^2');
end

function case_user_options()
% The options 'exclude' and 'degs' of the user are not overwritten
pvar s
[prog,x] = new_prog({'s'});
f = innerprod(x,x);
[~,~,~,Zop0] = piesos_ineq(prog,f);
assert(size(Zop0,1)==1,'Default: Zop must have 1 row.');
% Excluding the multiplier term: 2 rows, and int x^2 (delta term) cannot be represented
[progE,~,~,ZopE] = piesos_ineq(prog,f,struct('exclude',[true;false;false]));
assert(size(ZopE,1)==2,'exclude(1): Zop must have 2 rows, got %d.',size(ZopE,1));
assert(~feasible(progE),'Without multiplier term, int x^2 cannot be represented.');
% For the nonlocal kernel, including the multiplier term gives 3 rows
K = intop({polynomial(2)},{[0 3 0 1]},{'s'},[0,1]);
[~,~,~,ZopI] = piesos_ineq(prog,make_V(x,K,2),struct('exclude',[false;false;false]));
assert(size(ZopI,1)==3,'exclude=false: Zop must have 3 rows, got %d.',size(ZopI,1));
% degs=2 gives more rows than the default (f without delta term)
fN = make_V(x,K,2);
[~,~,~,ZopN0] = piesos_ineq(prog,fN);
[~,~,~,Zop2] = piesos_ineq(prog,fN,struct('degs',2));
assert(size(Zop2,1)>size(ZopN0,1),'degs=2 must give a larger basis than the default.');
end

function case_zero_basis()
% No monomial can contribute (odd degree, Newton polytope): f=0 is imposed,
% and g is the zero polynomial. int x^3 has coefficient 1, so infeasible
pvar s
[prog,x] = new_prog({'s'});
f3 = innerprod(x.*x,x);
assert(isequal(full(f3.degmat),3),'Test setup: expected monomial of degree 3.');
[prog,g,Qmat,Zop,Z] = piesos_ineq(prog,f3);
assert(isa(g,'polyopvar') && isempty(Qmat),'Unexpected outputs for empty basis.');
assert(~feasible(prog),'int x^3 = 0 must be infeasible.');
end

function case_constants()
% Degree-0 monomials: constants are stored as tensopvar in f, intop in g
pvar s
[~,x] = new_prog({'s'});
f2 = innerprod(x,x);
check_ineq(f2+2,x,true,'int x^2 + 2');
check_ineq(f2-2,x,false,'int x^2 - 2');
check_ineq(f2+innerprod(x.*x,x.*x)+1,x,true,'int x^2 + int x^4 + 1');
check_ineq(f2-innerprod(x.*x,x.*x)+1,x,false,'int x^2 - int x^4 + 1');
end

function case_psatz()
% Weights positive on [0,1] only: psatz multipliers are needed
pvar s
[~,x] = new_prog({'s'});
check_ineq(innerprod(x,s,x),x,true,'int s*x^2, psatz 1',struct('psatz',1));
check_ineq(innerprod(x,s*(1-s),x),x,true,'int s(1-s)*x^2, psatz 2',struct('psatz',2));
% Negative on part of the interval: infeasible regardless
check_ineq(innerprod(x,s-0.5,x),x,false,'int (s-1/2)*x^2, psatz 1',struct('psatz',1));
end

function case_mixed_dependence()
% States x1(s1) and x2(s1,s2): functionals in each, and their sum
pvar s1 s2
x1 = polyopvar('x1',polynomial({'s1'}),[0,1]);
x2 = polyopvar('x2',polynomial({'s1','s2'}),[0,1;0,1]);
prog = piesos_program([x1;x2]);
fa = innerprod(x1,x1);
fb = innerprod(x2,x2);
for sgn = [1,-1]
    expected = (sgn==1);
    [progk,g,~,~,Z] = piesos_ineq(prog,sgn*(fa+fb));
    assert(numel(Z.varname)==2,'Z must be defined in both states.');
    assert(feasible(progk)==expected,'Sum of squares of states in different variables: wrong feasibility (sign %d).',sgn);
end
% Indefinite combination
[progk,~] = piesos_ineq(prog,fa-fb);
assert(~feasible(progk),'int x1^2 - int x2^2 must be infeasible.');

end

function case_vector_states()
% States of size larger than 1, with a scalar state: the basis of g has the degree per state (dmonomials
% has the states as variables), so for degree 1 the vectors x1 (size 2) and x2, one block of Z per state.
% The weights with cross terms between the entries of x1 (and between the states) are certified by the
% Gram matrices of these vectors
% NOTE: f of degree 2 in a vector-valued state (monomials x_i*x_j of its rows, f quartic) fails in
% intop/combine_terms ("Vector dimension of the distributed monomial does not match ..."), not tested
x = polyopvar({'x1';'x2'},[2;1],{'s'},[0,1],true(2,1));         % x1 of size 2, x2 scalar
[~,~,Qmat,~,Zb] = piesos_ineq(piesos_program(x),innerprod(x,x));
assert(isequal(sort(reshape(Zb.dims.out,1,[])),[1,2]),'The basis must have a block per state (sizes 1 and 2), got %s.',mat2str(Zb.dims.out));
assert(isequal(sort(cellfun(@(q) size(q,1),Qmat(:,1))'),[1,2]),'Qmat must have blocks of the sizes of the states.');
check_ineq(innerprod(x,x),x,true,'vector state, int |x|^2');
check_ineq(-innerprod(x,x),x,false,'vector state, -int |x|^2');
Ppd = [2 -1 .5; -1 2 .3; .5 .3 1];                              % positive definite, with cross terms
check_ineq(innerprod(x,Ppd,x),x,true,'vector state, int x''*P*x, P>0');
check_ineq(-innerprod(x,Ppd,x),x,false,'vector state, -int x''*P*x, P>0');
Pind = [1 2 0; 2 1 0; 0 0 1];                                   % indefinite in the rows of x1
check_ineq(innerprod(x,Pind,x),x,false,'vector state, int x''*P*x, P indefinite');
Ps = [1 0 -.9; 0 1 0; -.9 0 1];                                 % couples a row of x1 and the scalar state
check_ineq(innerprod(x,Ps,x),x,true,'vector and scalar state, int x''*P*x, P>0');
check_ineq(innerprod(x,Ps-.05*eye(3),x),x,true,'vector and scalar state, P>0 (eigenvalues .05, 1, 1.9)');
check_ineq(innerprod(x,Ps-.2*eye(3),x),x,false,'vector and scalar state, P indefinite (eigenvalue -.1)');
end

function case_two_variables()
pvar s1 s2
[~,y] = new_prog({'s1','s2'},'y');
check_ineq(innerprod(y,y),y,true,'int y^2, 2 variables');
check_ineq(-innerprod(y,y),y,false,'-int y^2, 2 variables');
check_ineq(innerprod(y,s1,y),y,true,'int s1*y^2, psatz [1,0]',struct('psatz',[1,0]));
check_ineq(innerprod(y,s1-0.5,y),y,false,'int (s1-1/2)*y^2, psatz [1,0]',struct('psatz',[1,0]));
end

function case_two_states()
% Two states in the same variable: the basis for degree 4 contains the cross monomial y1*y2,
% and zero-coefficient monomials left by dmonomials must not end up in the basis
pvar s
y1 = polyopvar('y1',s,[0,1]);
y2 = polyopvar('y2',s,[0,1]);
prog = piesos_program([y1;y2]);
% int y1^4 + y2^4 - 1.5 y1^2 y2^2 >= 0 (Gram matrix [1 -.75;-.75 1] of (y1^2,y2^2) is positive
% definite), but int y1^4 + y2^4 - 3 y1^2 y2^2 < 0 for y1=y2
fy = innerprod(y1.*y1,y1.*y1)+innerprod(y2.*y2,y2.*y2);
cr = innerprod(y1.*y1,y2.*y2);
[progk,~,~,~,Z] = piesos_ineq(prog,fy-1.5*cr);
assert(isequal(sortrows(full(Z.Z)),[1 1;1 2;2 2]),'Basis must be [y1^2;y1*y2;y2^2], got %s.',mat2str(full(Z.Z)));
assert(feasible(progk),'int y1^4 + y2^4 - 1.5 y1^2 y2^2 must be feasible.');
[progk,~] = piesos_ineq(prog,fy-3*cr);
assert(~feasible(progk),'int y1^4 + y2^4 - 3 y1^2 y2^2 must be infeasible.');
end

function case_uncombined_f()
% f whose monomials are not combined or pruned must give the basis of the same functional
% written without redundancy, int x1^2 (the redundancy only enlarges the SDP, not f).
% (a) duplicate rows x1^2 with kernels 1, s^4, -s^4 (built by hand: plus and times combine
%     terms, but do not remove a monomial left with a zero operator)
% (b) int x1^2 + int x2^4 - int x2^4 from plus/minus: the x2^4 row is kept with a zero
%     operator (combine_terms does not prune it)
pvar s
x1 = polyopvar('x1',polynomial({'s'}),[0,1]);
x2 = polyopvar('x2',polynomial({'s'}),[0,1]);
prog = piesos_program([x1;x2]);
fA = innerprod(x1,x1);                  % int x1^2
fB = innerprod(x1,s^4,x1);              % int s^4 x1^2
[~,~,Qref,Zref,Zbas_ref] = piesos_ineq(prog,fA);
assert(isequal(block_sizes(Qref),1) && isequal(size(Zref),[1,1]),'Test setup: reference basis must be 1 x 1 (int x1^2, delta term only).');

% (a) Duplicate monomial rows
f = stack_rows({fA,fB,-fB});
assert(size(f.Z,1)==3 && size(combine_terms(f).Z,1)==1,'Test setup: f must have 3 rows, combining to 1.');
[progk,~,Qmat,Zop,Zbas] = piesos_ineq(prog,f);
assert(isequal(block_sizes(Qmat),block_sizes(Qref)) && isequal(size(Zop),size(Zref)), ...
    '(a) Basis for duplicated rows must be 3 x 3 as for int x1^2, got Qmat %s, Zop %s.',mat2str(block_sizes(Qmat)),mat2str(size(Zop)));
assert(isequal(full(Zbas.Z),full(Zbas_ref.Z)),'(a) Distributed monomial basis must be that of int x1^2.');
assert(feasible(progk),'(a) int x1^2 (written with cancelling duplicates) must be feasible.');
f = stack_rows({-fA,fB,-fB});
[progk,~] = piesos_ineq(prog,f);
assert(~feasible(progk),'(a) -int x1^2 (written with cancelling duplicates) must be infeasible.');

% (b) Monomial with zero operator
g4 = innerprod(x2.*x2,x2.*x2);          % int x2^4
for sgn = [1,-1]
    f = sgn*fA + g4 - g4;
    assert(size(f.Z,1)==2,'Test setup: f must have the monomials x1^2 and x2^4.');
    [progk,~,Qmat,Zop,Zbas] = piesos_ineq(prog,f);
    if sgn==1
        assert(isequal(block_sizes(Qmat),block_sizes(Qref)) && isequal(size(Zop),size(Zref)), ...
            '(b) Basis must be 3 x 3 as for int x1^2, got Qmat %s, Zop %s.',mat2str(block_sizes(Qmat)),mat2str(size(Zop)));
        assert(numel(Zbas.varname)==1,'(b) The basis must not involve x2.');
    end
    assert(feasible(progk)==(sgn==1),'(b) int x1^2 + (int x2^4 - int x2^4): wrong feasibility (sign %d).',sgn);
end
end

function f = stack_rows(fs)
% Scalar functional with the monomial row of each f in fs stacked, without combine_terms
f = fs{1};
ops = cellfun(@(g) g.C.ops{1,1},fs,'UniformOutput',false);
f.Z = sparse(cell2mat(cellfun(@(g) full(g.Z),fs(:),'UniformOutput',false)));
f.C = tensopmat(ops);
f.chkval = false;
f = initialize(f);
end

function sz = block_sizes(Qmat)
% Row dimensions of the blocks of Qmat
sz = cellfun(@(q) size(q,1),Qmat);
end

function case_mixed_basis_warning()
% A monomial in the basis for g cannot mix states depending on different variables
% (piesos_quadvar): piesos_ineq discards it and must warn, as f>=0 is then conservative.
% The cross term between the basis monomials x1 and x2 does not mix (only a product of two
% basis monomials does): no warning.
% (a) f = int int (x1*x2)^2 for x1(s1), x2(s1,s2): the only basis monomial x1*x2 mixes, so
%     f=0 is imposed (infeasible, although f is a square): the known limitation
% (b) f = int x1^2 + int int (x1*x2)^2: basis [x1; x1*x2], only the latter mixes
% (c) no warning: x1 and x2 in the same variable, or basis [x1;x2]
id = 'PIESOS:ineq:mixedDependence';
prev = warning('query',id);  warning('on',id);
restore = onCleanup(@() warning(prev));
x1 = polyopvar('x1',polynomial({'s1'}),[0,1]);
x2 = polyopvar('x2',polynomial({'s1','s2'}),[0,1;0,1]);
prog = piesos_program([x1;x2]);

% (a) Only monomial mixes
f = innerprod(x1.*x2,x1.*x2);
assert(isequal(full(f.degmat),[2,2]),'Test setup: f must be of degree [2,2] in [x1,x2].');
lastwarn('');
evalc('[progk,~,Qmat] = piesos_ineq(prog,f);');
[msg,wid] = lastwarn;
assert(strcmp(wid,id),'(a) Expected warning %s, got "%s".',id,wid);
assert(contains(msg,'1 of 1') && contains(msg,'f=0 is imposed'),'(a) Warning must state that no monomial remains: "%s".',msg);
assert(isempty(Qmat),'(a) Basis must be empty.');
assert(~feasible(progk),'(a) f=0 is imposed, so int int (x1*x2)^2 must be infeasible (known limitation).');

% (b) One of two monomials mixes
f = innerprod(x1,x1) + innerprod(x1.*x2,x1.*x2);
assert(isequal(sortrows(full(f.degmat)),[2,0;2,2]),'Test setup: f must have the monomials x1^2 and x1^2*x2^2.');
lastwarn('');
evalc('[~,~,Qmat,~,Z] = piesos_ineq(prog,f);');
[msg,wid] = lastwarn;
assert(strcmp(wid,id),'(b) Expected warning %s, got "%s".',id,wid);
assert(contains(msg,'1 of 2') && ~contains(msg,'f=0 is imposed'),'(b) Warning must state that 1 of 2 monomials is discarded: "%s".',msg);
assert(isequal(full(Z.Z),1) && numel(Z.varname)==1 && isequal(size(Qmat),[1,1]),'(b) Basis must be [x1].');

% (c) No warning for monomials that do not mix
y1 = polyopvar('y1',polynomial({'s'}),[0,1]);
y2 = polyopvar('y2',polynomial({'s'}),[0,1]);
progy = piesos_program([y1;y2]);
f = innerprod(y1.*y2,y1.*y2);           % same dependence: basis [y1*y2] does not mix
lastwarn('');
evalc('[progk,~,~,~,Z] = piesos_ineq(progy,f);');
[~,wid] = lastwarn;
assert(~strcmp(wid,id),'(c) No warning expected for monomials of states in the same variable.');
assert(isequal(full(Z.Z),[1,2]) && feasible(progk),'(c) int (y1*y2)^2 must be feasible with basis [y1*y2].');
lastwarn('');
evalc('[~,~,~,~,Z] = piesos_ineq(prog,innerprod(x1,x1)+innerprod(x2,x2));');
[~,wid] = lastwarn;
assert(~strcmp(wid,id),'(c) No warning expected for the basis [x1;x2].');
assert(numel(Z.varname)==2,'(c) Basis must involve both states.');
end

function case_kernel_degrees()
% Degrees d (dummy variables) and p (primary variable) of Zop, and terms in Zop, from the
% exponents of the kernels of f, one variable s on [0,1]. For the kernel monomial s1^a*s2^b
% of a quadratic f, g can only have the monomials with min(a,b)<=d and max(a,b)<=d+2p+1, so
% d and p are needed and sufficient. Oracle: (a) f = (int s^a x)^2; (b) f SOS by
% construction, int |W(s)*Zop*x|^2, Zop of known degree dd0 and W of degree dp0 (monopmat,
% innerprod, random W), for which f>=0 holds for any (dd0,dp0)
pvar s s_dum1 s_dum2
[~,x] = new_prog({'s'});
prog = piesos_program(x);
% (a) d=a is needed, not the total degree 2a of the kernel
for a=0:3
    K = intop({2*s_dum1^a*s_dum2^a},{[0 3 0 1]},{'s'},[0,1]);
    f = make_V(x,K,2);
    [~,~,~,Zop] = piesos_ineq(prog,f);
    assert(size(Zop,1)==2*(a+1),'(int s^%d x)^2: Zop must have %d rows (no multiplier), got %d.',a,2*(a+1),size(Zop,1));
    check_ineq(f,x,true,sprintf('(int s^%d x)^2',a));
    if a>0
        check_ineq(f,x,false,sprintf('(int s^%d x)^2, degs=%d',a,a-1),struct('degs',a-1));
    end
end
% (b) SOS by construction: feasible by default with Zop of the known size, infeasible for
% smaller d or p
for c = [0 0; 0 2; 1 1; 2 1; 3 0]'
    dd0 = c(1);     dp0 = c(2);
    label = sprintf('SOS, dd0=%d, dp0=%d',dd0,dp0);
    [f,xs] = make_sos_f1(dd0,dp0,[0;0;0],5);
    progs = piesos_program(xs);
    [~,~,~,Zop] = piesos_ineq(progs,f);
    assert(size(Zop,1)==2*dd0+3,'%s: Zop must have %d rows, got %d.',label,2*dd0+3,size(Zop,1));
    check_ineq(f,xs,true,label);
    if dd0>0
        check_ineq(f,xs,false,[label,', d-1'],struct('degs',struct('in',dd0-1,'out',dp0)));
    end
    if dp0>0
        check_ineq(f,xs,false,[label,', p-1'],struct('degs',struct('in',dd0,'out',dp0-1)));
    end
end
% (c) Delta term only: Zop = x, p from the degree of the weight: int s^2*x^2
% is the square int (s*x)^2, needs p=1 and no integral terms
f = innerprod(x,s^2,x);
[~,~,Q,Zop] = piesos_ineq(prog,f);
assert(size(Zop,1)==1 && isequal(size(Q{1,1}),[1,1]),'int s^2 x^2: Zop must be x (1 row).');
check_ineq(f,x,true,'int s^2 x^2');
check_ineq(f,x,false,'int s^2 x^2, p=0',struct('degs',struct('in',0,'out',0)));
% (d) Delta term and nonlocal kernel (1+s*t)^2 = 1 + 2st + s^2t^2: d=2, all terms
K = intop({2*(1+s_dum1*s_dum2)^2},{[0 3 0 1]},{'s'},[0,1]);
f = innerprod(x,x)+make_V(x,K,2);
[~,~,~,Zop] = piesos_ineq(prog,f);
assert(size(Zop,1)==7,'int x^2 + int int (1+st)^2 xx: Zop must have 7 rows, got %d.',size(Zop,1));
check_ineq(f,x,true,'int x^2 + int int (1+st)^2 x(s)x(t)');
% (e) Degree 4, delta terms only: Zop = x^2 (1 row), and for x^2 + x^4 [x; x^2]
f = innerprod(x.*x,x.*x);
[~,~,~,Zop] = piesos_ineq(prog,f);
assert(size(Zop,1)==1,'int x^4: Zop must have 1 row, got %d.',size(Zop,1));
check_ineq(f,x,true,'int x^4');
f = innerprod(x,x)+innerprod(x.*x,x.*x);
[~,~,~,Zop] = piesos_ineq(prog,f);
assert(size(Zop,1)==2,'int x^2 + int x^4: Zop must have 2 rows, got %d.',size(Zop,1));
check_ineq(f,x,true,'int x^2 + int x^4');
end

function case_option_types()
% Degrees and terms per spatial variable, two variables s1,s2 on [0,1]: f SOS by construction
% with the term types (1: multiplier, 2: lower, 3: upper integral) allowed in each variable
% and degrees dd0, dp0 per variable. piesos_ineq must infer these: Zop of the known size,
% feasible, and infeasible for a smaller degree
ALL = [1 2 3];
cases = {[0 0],[0 0],{1,1};         % delta terms in both variables
         [0 0],[0 0],{1,[2 3]};     % local in s1, nonlocal in s2
         [1 1],[0 0],{[2 3],1};     % nonlocal in s1, local in s2
         [1 0],[0 0],{[2 3],[2 3]}; % no multiplier in either variable
         [1 0],[0 1],{ALL,ALL};     % all terms
         [0 0],[1 0],{ALL,ALL}};
rows = @(d,al) prod(cellfun(@(a,dv) nnz(a==1)+nnz(a>1)*(dv+1),al,num2cell(d)));
for c=1:size(cases,1)
    dd0 = cases{c,1};   dp0 = cases{c,2};   al = cases{c,3};
    label = sprintf('d=%s, p=%s, types %s|%s',mat2str(dd0),mat2str(dp0),mat2str(al{1}),mat2str(al{2}));
    [f,y] = make_sos_f2(dd0,dp0,al,6);
    prog = piesos_program(y);
    [~,~,~,Zop] = piesos_ineq(prog,f);
    % rows: (1 if multiplier + 2(d+1) integral terms) per variable; f determines the degree
    % only of variables with an integral term
    nr = rows(dd0,al);
    assert(size(Zop,1)==nr,'%s: Zop must have %d rows, got %d.',label,nr,size(Zop,1));
    check_ineq(f,y,true,label);
    for v=1:2
        if dd0(v)>0 && any(al{v}>1)
            dd = dd0;   dd(v) = dd(v)-1;
            check_ineq(f,y,false,sprintf('%s, d(%d)-1',label,v),struct('degs',struct('in',dd,'out',dp0)));
        end
        if dp0(v)>0
            pp = dp0;   pp(v) = pp(v)-1;
            check_ineq(f,y,false,sprintf('%s, p(%d)-1',label,v),struct('degs',struct('in',dd0,'out',pp)));
        end
    end
end
end

function case_psatz_degree()
% With psatz multipliers of odd order q, Pmat has degree q+2(p-ceil(q/2)): a weight of degree 2
% on [0,1] needs p=2 for q=1 (s*sigma1 + (1-s)*sigma2, leading terms cancelling), not p=1:
% int ((s-0.3)^2+0.1) x^2, psatz 1
pvar s
[~,x] = new_prog({'s'});
f = innerprod(x,(s-0.3)^2+0.1,x);
check_ineq(f,x,true,'int ((s-.3)^2+.1) x^2, psatz 1',struct('psatz',1));
check_ineq(f,x,false,'int ((s-.3)^2+.1) x^2, psatz 1, p=1',struct('psatz',1,'degs',struct('in',0,'out',1)));
check_ineq(f,x,true,'int ((s-.3)^2+.1) x^2, psatz 1, p=2',struct('psatz',1,'degs',struct('in',0,'out',2)));
% Without psatz, p=1 suffices
check_ineq(f,x,true,'int ((s-.3)^2+.1) x^2, p=1',struct('degs',struct('in',0,'out',1)));
end

function case_symbolic_f()
% f with decision variables: V (SOS, d=p=1) + int x^2. Monomials with only decision variables
% need not be recovered, but d=p=ceil((e(1)-1)/3)=1 is used for them (e(1)=4)
pvar s
[~,x] = new_prog({'s'});
prog = piesos_program(x);
[prog,V] = piesos_sosvar(prog,dmonomials(x,1),struct('degs',struct('in',1,'out',1)));
f = V + innerprod(x,x);
[progk,~,~,Zop] = piesos_ineq(prog,f);
assert(size(Zop,1)==5,'f with decision variables: Zop must have 5 rows (d=1, all terms), got %d.',size(Zop,1));
assert(feasible(progk),'V + int x^2 with V SOS must be feasible.');
end

function [f,x] = make_sos_f1(dd0,dp0,excl,seed)
% f(x) = int_0^1 |W(s)*Zop*x(s)|^2 ds, Zop*x = [x; int_0^s z*x; int_s^1 z*x]
% (types not excluded by excl) with z = [1,...,s^dd0], and W a random matrix
% polynomial of degree dp0: f>=0 by construction, with the SOS of the given degrees
pvar s
rng(seed);
x = polyopvar('x1',s,[0,1]);
Zop = monopmat(1,struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0 1],'out',[0 1]), ...
    struct('in',1,'out',true),struct('in',dd0,'out',0),struct('exclude',excl(:)));
r = size(Zop,1);
W = polynomial(zeros(r,r));
for k=0:dp0
    W = W + s^k*randn(r,r);
end
f = innerprod(Zop*x,W'*W,[]);
end

function [f,y] = make_sos_f2(dd0,dp0,allowed,seed)
% As make_sos_f1, for y(s1,s2): the terms of Zop are the products of the types
% allowed{v} in variable v (1: multiplier, 2: lower, 3: upper integral), with
% dummy degree dd0(v), and W is of degree dp0(v) in s_v
pvar s1 s2
rng(seed);
y = polyopvar('y1',polynomial({'s1','s2'}),[0,1;0,1]);
ex = true(3,3);
ex(allowed{1},allowed{2}) = false;
Zop = monopmat(1,struct('in',{{'s1','s2'}},'out',{{'s1','s2'}}),struct('in',[0 1;0 1],'out',[0 1;0 1]), ...
    struct('in',[1 1],'out',[true true]),struct('in',dd0(:)','out',[0 0]),struct('exclude',ex));
r = size(Zop,1);
W = polynomial(zeros(r,r));
for k1=0:dp0(1)
    for k2=0:dp0(2)
        W = W + s1^k1*s2^k2*randn(r,r);
    end
end
f = innerprod(Zop*y,W'*W,[]);
end
