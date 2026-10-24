function Test_piesos_eq()
% TEST_PIESOS_EQ tests, using assertions, piesos_eq,
%   prog = piesos_eq(prog,f),
% which adds the constraint f=0 to a PIESOS program for a scalar
% 'polyopvar' f whose coefficients are 'intop' (or 'dpvar') objects.
%
% Oracle, independent of the function under test:
%  - number of constraints: the monomial x o x is symmetric, so the
%    integrals over s_dum1<=s_dum2 and s_dum1>=s_dum2 act as one. 'full'
%    kernels in N variables, with 3^N terms (multiplier, lower, upper
%    along each variable), are equivalent to (3^N+1)/2 terms: the swap of
%    the two factors maps lower <-> upper along every variable, fixing
%    only the all-multiplier term. One constraint is added per term left
%    after combine_terms;
%  - values: Kop*(x o x) = T(x o x), for a target functional T of known
%    polynomial kernel kappa, is feasible exactly if the kernels can
%    represent T, and the solved functional is compared with the exact
%    integral of kappa*p(s)*p(s), from the moments
%    int_0^1 s^n ds = 1/(n+1) (not with apply_op); an unrepresentable
%    target must be infeasible;
%  - dpvar coefficients: the constraint (c1-3)=0 must give c1=3;
%  - workflow: V from piesos_sosvar set equal to innerprod(x,x), solved,
%    with the quadratic form of the solved Pmat integrated by quadrature.
%  - equal functionals written in different forms (case_canonical): the
%    constraint must not over-constrain the coefficients.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_eq
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
% DJ (CC), 10/23/2026: Merge Test_piesos_eq_canonical into this test, as
%   cases;

% NOTES
% - For N>=2 combine_terms does not merge all swapped terms (2x4 terms
%   for N=2, 2x4x4 for N=3, the same as for the functional before
%   merging), so only #constraints >= (3^N+1)/2 is asserted; N=1 gives
%   exactly (3^N+1)/2;
% - The kernels are not determined by the constraints (kappa(s1,s2) and
%   kappa(s2,s1) only enter through their sum), so only the solved
%   functional is compared;
% - Without a target the kernels solve to zero, which alone cannot
%   distinguish piesos_eq from a function adding no constraints;
%   case_intop_zero also asserts the constraint count;
% - Equality-only programs fail in SeDuMi ("Size b mismatch"): solve_ok
%   adds a trivial PSD matrix.

tests ={@case_input_errors, @case_dpvar_coefficient, @case_lumped_intop, ...
         @case_intop_zero, @case_intop_target_1var, @case_intop_infeasible, ...
         @case_intop_target_2var, @case_term_counts, @case_vartable, ...
         @case_zero_kernel, @case_sosvar_pipeline, ...
         @case_canonical};
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
    error('Test_piesos_eq: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_eq: all %d cases passed.\n',numel(tests));
end

%% ---- cases ------------------------------------------------------------
function case_canonical()
% case_canonical -- tests that piesos_eq does not
% over-constrain the coefficients of equal functionals. The equation
%   Kop*(x o x) = T*(x o x),
% for T a given functional on the full square [0,3,0,3] with a kernel that
% is not symmetric, and Kop the unknown functional of piesos_intopvar with
% all integrals (multiplier, lower [0,3,0,1], and upper [0,2,0,3]), must
% be feasible, and the solved functional must equal T.
%
% The unknown functional has no integral over the full square, only the
% two chains: the equation is feasible only if combine_terms writes the
% square as the chain [0,3,0,1] with the kernel k(s1,s2)+k(s2,s1). This is
% also tested for T given by the lower chain and the kernel
% k(s1,s2)+k(s2,s1), and by the upper chain, which give the same T.
%

% NOTES
% - Only the solved functional is compared with T, on random products
%   p(s_1)*p(s_2): the kernels of Kop are not unique.

    pvar s s_dum1 s_dum2
    rng(3);
    kap = 1+2*s_dum1-s_dum2+0.5*s_dum1*s_dum2;          % not symmetric
    ksw = subs(kap,[s_dum1;s_dum2],[s_dum2;s_dum1]);
    kmul = 2-s_dum1+s_dum1^2;

    % T given in three ways: the square, the lower chain, and the upper chain
    Tarr = cell(1,3);
    Tarr{1} = intop({kmul;kap},{[0,3,1,1;0,3,0,3]},{'s'},[0,1]);
    Tarr{2} = intop({kmul;kap+ksw},{[0,3,1,1;0,3,0,1]},{'s'},[0,1]);
    Tarr{3} = intop({kmul;kap+ksw},{[0,3,1,1;0,2,0,3]},{'s'},[0,1]);

    for j=1:3
        x1 = polyopvar('x1',polynomial({'s'}),[0,1]);
        prog = piesos_program(x1);
        [prog,Kop] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',2);
        prog = piesos_eq(prog,make_V(x1,Kop-Tarr{j},2));
        [prog,~] = sosposmatr(prog,1);      % an equality-only program cannot be solved by SeDuMi
        out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
        if ~(info.pinf==0 && info.dinf==0)
            error("case_canonical failed: form %d of T, the equation Kop = T is infeasible.",j)
        end
        Ks = piesos_getsol(prog,Kop);
        for t=1:3
            c = randn(3,1);
            pvar s_1 s_2
            xval = (c(1)+c(2)*s_1+c(3)*s_1^2)*(c(1)+c(2)*s_2+c(3)*s_2^2);
            val = double(apply_op(Ks,xval));
            ref = double(apply_op(Tarr{j},xval));
            if abs(val-ref)>1e-6*max(1,abs(ref))
                error("case_canonical failed: form %d of T, the solved functional differs from T: %g vs %g.",j,val,ref)
            end
        end
    end

    disp('case_canonical passed.')
end

function [prog,x1] = new_prog(varnames)
% Program in the distributed state x1 of the given spatial variables, on [0,1]
x1 = polyopvar('x1',polynomial(varnames(:)),[0,1]);
prog = piesos_program(x1);
end

function V = make_V(x1,K,d)
% d-th power of the state, with coefficient functional K, as in piesos_sosvar
V = x1;
for i=2:d
    V = V*x1;
end
V.C = tensopmat({K});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
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

function [is_ok,prog] = solve_ok(prog)
[prog,~] = sosposmatr(prog,1);
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
is_ok = (info.pinf==0 && info.dinf==0);
end

function names = dummy_names(vnames,deps)
% Names s_dumj of the dummy variables of the kernels, ordered as the oracle exponents:
% all copies of the first variable, then of the second, etc.
names = cell(1,0);
for i=1:numel(vnames)
    for j=1:deps
        names{end+1} = [vnames{i},'_dum',num2str(j)];      %#ok<AGROW>
    end
end
end

function kap = kernel_poly(terms,dnames)
% 'polynomial' sum_t terms(t,1)*prod_k dnames{k}^terms(t,1+k)
sv = polynomial(dnames(:));
kap = polynomial(0);
for t=1:size(terms,1)
    mon = polynomial(terms(t,1));
    for k=1:numel(dnames)
        mon = mon*sv(k)^terms(t,1+k);
    end
    kap = kap+mon;
end
end

function p = poly_from_array(P,vnames)
% 'polynomial' sum_k P(k)*prod_i vnames{i}^(j_i(k)), for the N-d array P,
% with (j_1,...,j_N) the zero-based position of P(k)
N = numel(vnames);
sz = [size(P),ones(1,N-ndims(P))];
sv = polynomial(vnames(:));
p = polynomial(0);
idx = cell(1,N);
for k=1:numel(P)
    [idx{:}] = ind2sub(sz,k);
    mon = polynomial(P(k));
    for i=1:N
        mon = mon*sv(i)^(idx{i}-1);
    end
    p = p+mon;
end
end

function val = moment_oracle(terms,P)
% Exact value of int_{[0,1]^(2N)} kappa(s_dum)*p(copy 1)*p(copy 2) over all dummy variables,
% for kappa given by terms (see kernel_poly, dummy variables as in dummy_names with 2 copies)
% and p by P (see poly_from_array), from the moments int_0^1 s^n ds = 1/(n+1)
N = (size(terms,2)-1)/2;
sz = [size(P),ones(1,N-ndims(P))];
J = zeros(numel(P),N);
idx = cell(1,N);
for k=1:numel(P)
    [idx{:}] = ind2sub(sz,k);
    J(k,:) = [idx{:}]-1;
end
val = 0;
for t=1:size(terms,1)
    e = terms(t,2:end);
    for k1=1:numel(P)
        for k2=1:numel(P)
            w = terms(t,1)*P(k1)*P(k2);
            for i=1:N
                w = w/(e(2*i-1)+J(k1,i)+1)/(e(2*i)+J(k2,i)+1);
            end
            val = val+w;
        end
    end
end
end

function val = diag_oracle(mterms,P)
% Exact value of int_0^1 kmul(t)*p(t)^2 dt for the polynomial p of one
% variable with coefficients P(j+1) of t^j, and kmul = sum_t mterms(t,1)*t^mterms(t,2)
val = 0;
for t=1:size(mterms,1)
    for j=1:numel(P)
        for k=1:numel(P)
            val = val+mterms(t,1)*P(j)*P(k)/(mterms(t,2)+j+k-1);
        end
    end
end
end

function case_input_errors()
[prog,x1] = new_prog({'s'});
assert_err(@() piesos_eq(prog,3),'polyopvar');
% coefficients of x1*x1 are 'tensopvar' (the identity): the call must fail on these
assert_err(@() piesos_eq(prog,x1*x1),'should be specified as ''dpvar'' or ''intop''');
end

function case_dpvar_coefficient()
% 'dpvar' coefficients, assigned to the property as a 'tensopmat' cannot declare them:
% the constraint c1-3=0 must be added, and give c1=3
[prog,x1] = new_prog({'s'});
dpvar c1;
prog = sosdecvar(prog,c1);
V = x1*x1;
C = V.C;
C.ops = {c1-3};
V.C = C;
V.chkval = true;
n0 = prog.expr.num;
prog = piesos_eq(prog,V);
assert(prog.expr.num==n0+1,'Expected one constraint, got %d.',prog.expr.num-n0);
[ok,prog] = solve_ok(prog);
assert(ok,'Program must be feasible.');
assert(abs(double(piesos_getsol(prog,c1))-3)<1e-6,'Constraint must give c1=3.');
end

function case_lumped_intop()
% State without spatial variables: constant 'intop' coefficients with decision variables,
% for 2 monomials y o y and y o y o y. Each monomial must be constrained
y = polyopvar('y');
prog = piesos_program(y);
dpvar c1 c2;
prog = sosdecvar(prog,[c1;c2]);
V = y*y + y*y*y;
assert(isequal(full(V.Z),[1 1 0;1 1 1]),'Unexpected order of the monomials.');
V.C = tensopmat({intop(c1-3),intop(c2+1)});
V.dims.out = V.C.dims.out;      V.vars.out = V.C.vars.out;
V.dom.out = V.C.dom.out;        V.depmat.out = V.C.depmat.out;
V.chkval = true;
n0 = prog.expr.num;
prog = piesos_eq(prog,V);
assert(prog.expr.num==n0+2,'Expected two constraints, got %d.',prog.expr.num-n0);
[ok,prog] = solve_ok(prog);
assert(ok,'Program must be feasible.');
assert(abs(double(piesos_getsol(prog,c1))-3)<1e-6 && abs(double(piesos_getsol(prog,c2))+1)<1e-6,...
    'Constraints must give c1=3, c2=-1.');
end

function case_intop_zero()
% The functional of the 'full' kernels of degree 2 and cumulative degree 1,
% acting on x o x, set to zero: 2 constraints, and the kernels solve to zero
[prog,x1] = new_prog({'s'});
[prog,Kop] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',1);
V = make_V(x1,Kop,2);
n0 = prog.expr.num;
prog = piesos_eq(prog,V);
assert(prog.expr.num==n0+2,'Expected 2 constraints, got %d.',prog.expr.num-n0);
[ok,prog] = solve_ok(prog);
assert(ok,'Program must be feasible.');
Ks = piesos_getsol(prog,Kop);
Ks_kernels = kernels(Ks);
for k=1:numel(Ks.params)
    assert(isempty(Ks_kernels{k}.C) || max(abs(Ks_kernels{k}.C(:)))<1e-8,'Kernel %d must be zero.',k);
end
end

function case_intop_target_1var()
% Kop*(x o x) = T(x o x), T the functional with kernel kappa on [0,1]^2, and
% kernel kmul on the diagonal s_dum1=s_dum2: the solved functional must equal
% T for all x, which constrains the multiplier and the integral term
rng(1);
terms = [1 0 0; 2 1 0; -1 0 1; 0.5 1 1];            % 1+2*s_dum1-s_dum2+0.5*s_dum1*s_dum2
mterms = [2 0; -1 1; 1 2];                          % 2-s_dum1+s_dum1^2
[prog,x1] = new_prog({'s'});
[prog,Kop] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',2);
assert(isequal(Kop.intmat{1},[0,3,1,1; 0,3,0,1; 0,2,0,3]),'Order of the terms: multiplier, lower, upper expected.');
kap = kernel_poly(terms,dummy_names({'s'},2));
kmul = kernel_poly(mterms,{'s_dum1'});
Kt = intop({kmul; kap; kap},Kop.intmat,Kop.vars,Kop.dom);   % kappa on s_dum2<=s_dum1 and s_dum1<=s_dum2
prog = piesos_eq(prog,make_V(x1,Kop-Kt,2));
[ok,prog] = solve_ok(prog);
assert(ok,'The functional with kernel kappa must be feasible.');
Ks = piesos_getsol(prog,Kop);
for t=1:4
    P = randn(4,1);                                 % random cubic p
    val = double(apply_op(Ks,poly_from_array(P,{'s_1'})*poly_from_array(P,{'s_2'})));
    ref = moment_oracle(terms,P)+diag_oracle(mterms,P);
    assert(abs(val-ref)<=1e-6*max(1,abs(ref)),'Solved functional differs from the target: %g vs %g.',val,ref);
end
end

function case_intop_infeasible()
% A target of degree 3 cannot be represented by kernels of degree 2
terms = [1 0 0; 1 3 0];
[prog,x1] = new_prog({'s'});
[prog,Kop] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',2);
kap = kernel_poly(terms,dummy_names({'s'},2));
Kt = intop({polynomial(0); kap; kap},Kop.intmat,Kop.vars,Kop.dom);
prog = piesos_eq(prog,make_V(x1,Kop-Kt,2));
assert(~solve_ok(prog),'A kernel of degree 3 cannot be represented with degree 2.');
end

function case_intop_target_2var()
% Two spatial variables, 3x3 kernels (multiplier, lower, upper along each):
% at least 5 constraints, and the solved functional equals the target
rng(2);
% kappa(a,b,c,d), (a,b)=(s1_dum1,s1_dum2), (c,d)=(s2_dum1,s2_dum2):
% 1+2a-b-d+0.5*a*c+b*d, of degree 1 in each of the two directions
terms = [1 0 0 0 0; 2 1 0 0 0; -1 0 1 0 0; -1 0 0 0 1; 0.5 1 0 1 0; 1 0 1 0 1];
vn = {'s1','s2'};
[prog,x1] = new_prog(vn);
[prog,Kop] = piesos_intopvar(prog,[1,1],vn,[0 1;0 1],[2,2],'full',{1,1});
kap = kernel_poly(terms,dummy_names(vn,2));
Kt_kernels = repmat({polynomial(0)},[3,3]);
Kt_kernels(2:3,2:3) = {kap};                         % lower/upper along both variables
Kt = intop(Kt_kernels,Kop.intmat,Kop.vars,Kop.dom);
n0 = prog.expr.num;
prog = piesos_eq(prog,make_V(x1,Kop-Kt,2));
ncon = prog.expr.num-n0;
assert(ncon>=5,'Expected at least 5 constraints, got %d.',ncon);
[ok,prog] = solve_ok(prog);
assert(ok,'The functional with kernel kappa must be feasible.');
Ks = piesos_getsol(prog,Kop);
for t=1:3
    P = randn(3,3);                                 % random p(s1,s2) of degree 2 in each
    pv = poly_from_array(P,{'s1_1','s2_1'})*poly_from_array(P,{'s1_2','s2_2'});
    val = double(apply_op(Ks,pv));
    ref = moment_oracle(terms,P);
    assert(abs(val-ref)<=1e-6*max(1,abs(ref)),'Solved functional differs from the target: %g vs %g.',val,ref);
end
end

function case_term_counts()
% At least (3^N+1)/2 constraints for the 'full' kernels in N variables acting
% on x o x, and exactly (3^N+1)/2 for N=1
for N=1:3
    vn = arrayfun(@(i) sprintf('s%d',i),1:N,'UniformOutput',false);
    [prog,x1] = new_prog(vn);
    [prog,Kop] = piesos_intopvar(prog,[1,1],vn,repmat([0 1],N,1),2*ones(1,N),'full',0);
    assert(numel(Kop.params)==3^N,'Expected 3^N terms.');
    n0 = prog.expr.num;
    prog = piesos_eq(prog,make_V(x1,Kop,2));
    ncon = prog.expr.num-n0;
    nmin = (3^N+1)/2;
    assert(ncon>=nmin,'N=%d: expected at least %d constraints, got %d.',N,nmin,ncon);
    assert(N>1 || ncon==nmin,'N=1: expected %d constraints, got %d.',nmin,ncon);
end
end

function case_vartable()
% The dummy variables of the kernels are added to the program once each, as a column
vn = {'s1','s2'};
[prog,x1] = new_prog(vn);
[prog,Kop] = piesos_intopvar(prog,[1,1],vn,[0 1;0 1],[2,2],'full',0);
V = make_V(x1,Kop,2);
prog = piesos_eq(prog,V);
expected = sort([vn,dummy_names(vn,2)]);
assert(iscolumn(prog.vartable) && isequal(sort(prog.vartable(:))',expected),...
    'Table of variables: %s.',strjoin(prog.vartable(:)',','));
% a second call adds no variables
prog2 = piesos_eq(prog,V);
assert(isequal(prog2.vartable,prog.vartable),'Variables must not be added twice.');
end

function case_sosvar_pipeline()
% SOS functional V = <Zop*x, Pmat*Zop*x> from piesos_sosvar (Zop*x = [x; int_0^s x; int_s^1 x]),
% set equal to innerprod(x,x). The solved Pmat must be PSD with Pmat(1,1)=1 (other entries are
% not unique: the last two components add up to a constant), and the quadratic form integrated
% against random p must give int p^2, by quadrature
rng(5);
pvar s;
x = polyopvar('x1',s,[0,1]);
prog = piesos_program(x);
[prog,V,Pmat] = piesos_sosvar(prog,dmonomials(x,1),struct('degs',0));
% Pmat: one block per monomial of Z, here 1 x 1 holding the 3 x 3 quadratic form
assert(iscell(Pmat) && isscalar(Pmat),'Pmat must be a 1 x 1 cell.');
Pmat = Pmat{1,1};
assert(isequal(size(Pmat),[3,3]),'Unexpected size of Pmat.');
prog = piesos_eq(prog,V-innerprod(x,x));
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
assert(info.pinf==0 && info.dinf==0,'Program must be feasible.');
P = double(piesos_getsol(prog,Pmat));
assert(min(eig((P+P')/2))>-1e-6,'Pmat must be positive semidefinite.');
assert(abs(P(1,1)-1)<1e-5,'Pmat(1,1) must be 1, got %g.',P(1,1));
for t=1:3
    c = randn(1,3);                                 % p(s) = c1+c2*s+c3*s^2
    p = @(t) c(1)+c(2)*t+c(3)*t.^2;
    Ip = @(t) c(1)*t+c(2)*t.^2/2+c(3)*t.^3/3;       % int_0^t p
    z = @(t) [p(t); Ip(t); Ip(1)-Ip(t)];
    val = integral(@(t) arrayfun(@(ti) z(ti)'*P*z(ti),t),0,1,'AbsTol',1e-12,'RelTol',1e-10);
    ref = integral(@(t) p(t).^2,0,1,'AbsTol',1e-12,'RelTol',1e-10);
    assert(abs(val-ref)<1e-5*max(1,abs(ref)),'V(p) must equal int p^2: %g vs %g.',val,ref);
end
end

function case_zero_kernel()
% A kernel that is zero, without decision variables (the multiplier term has no monomials, as
% every explicit monomial depends on the variable set equal to another one), is accepted
[prog,x1] = new_prog({'s'});
% Monomials given as struct field 'full': a 1 x N cell with the monomials of each variable
[prog,Kop] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',struct('full',{{[0 1; 0 2]}}));
Kop_kernels = kernels(Kop);
assert(~isa(Kop_kernels{1},'dpvar'),'Multiplier kernel must be a zero polynomial.');
prog = piesos_eq(prog,make_V(x1,Kop,2));
[ok,~] = solve_ok(prog);
assert(ok,'Program must be feasible.');
end
