function Test_PIESOS_edge_cases(only)
% TEST_PIESOS_EDGE_CASES(ONLY) tests edge cases of PIESOS found in the
% release review of 11/2026, each against an independent oracle:
%   case_eq_multirow        piesos_eq with a 2-row polynomial constrains both
%                           rows: (1-c1)*<x,x> = 0, (2-c2)*<x,x> = 0 gives
%                           c1 = 1, c2 = 2;
%   case_liediff_vector     Liediff of V = <u,u> for a vector-valued PDE
%                           state equals 2*<T*x,f(x)>;
%   case_intop_diff         d^0/ds^0 of a functional is the functional, d/ds
%                           is zero, for D = [0;1];
%   case_scalar_products    2*u'*f == 2*(u'*f), x*3 == 3*x; for a decision
%                           variable c, c times a functional and the
%                           weight form y'*c*y are accepted, c*x, x*c and
%                           c*xo (ODE state) give a clear error;
%   case_constructor        polyopvar refuses duplicate state names and
%                           domains [a,b] with a >= b;
%   case_otimes_params      otimes_params accepts a 'double' first and a
%                           cell-coefficient struct second, with the same
%                           product as in the other order;
%   case_mtimes_unchecked   mtimes_tensopvar accepts a valid 'intop' whose
%                           chkval is false;
%   case_quadpoly_poly      a 'quadPoly' times a 'polynomial' equals the
%                           polynomial product;
%   case_intop_plus         intop + nonzero constant gives a clear error;
%   case_exec_linear        the PIESOS executives accept a linear PIE
%                           ('pie_struct'): u_t = u_ss + r*u, u(0)=u(1)=0, is
%                           stable iff r < pi^2 (r = 5 verified, 10.5 not);
%   case_quadpoly_fixes     'quadPoly' concatenation with numeric inputs,
%                           integration over a variable the polynomial
%                           does not depend on, and linear indexing with a
%                           vector, against the same operations on
%                           'polynomial' objects;
%   case_dpvar_constant     V - c with a decision variable c: sums with
%                           numeric constants, the SOS bound
%                           max c s.t. int(u^4-u^2) - c >= 0, which is
%                           -1/4, and piesos_eq(V-c = V-2), c = 2;
%   case_known_function     x + s and s*x (s a spatial variable), against
%                           the integrals int(x*(x+s)), int(s*x^2) for
%                           x = 1+s^2, computed with 'polynomial';
%   case_unknown_options    a misspelled option gives one warning, valid
%                           options none;
%   case_minor_fixes        A.*T for a matrix A (with a single row of T
%                           repeated), kron(2,T), K(i,j).params of an
%                           'intop', piesos_eq(V-1) infeasible, and
%                           piesos_getsol of c*(1-V), and the largest c
%                           with 2*V + c*(1-V) >= 0 (c = 2), against
%                           values computed by hand;
%   case_functional_product x*(x'*x), (x'*x)*x, (1+x'*x)*x and
%                           otimes(x,x'*x) give the error naming a
%                           functional; xo*(x'*x) for an ODE state xo
%                           still equals 2*int(1+s^2)^2 at xo = 2;
%   case_vector_scalar      a vector times a scalar distributed polynomial
%                           ([u;u]*u, u*[u;u], [u;xo]*u, y*u for a vector
%                           state y, [y;u]*u, [V;V]*V for V = u'*u),
%                           against values computed by hand; y.*[u;u]
%                           (different row-blocks) gives a clear error;
%   case_quadopvar_order    piesos_quadopvar with vars = {'s2','s1'}:
%                           opts.intmat = {'steps','full'} follows this
%                           order as psatz does (steps in s2), and the
%                           returned intmat is in the order (s1,s2) of Zop;
%   case_apply_poly_checks  apply_poly with a tensor-product output errors
%                           for a missing state value, and for a value with
%                           decision variables;
%   case_review_operators   T*P for a 1x1 'sopvar' P composes as T(P(x)),
%                           [1,2].*T and [1;2].*T for a 2x2 T, V.*[V;V*V],
%                           [1;2]*(x.*x), T{:,1}, and an 'intop' times a
%                           'tensopmat' with an empty block;
%   case_review_polyopvar   subs of a constant into a functional of degree
%                           3, apply_poly with 'dpvar' coefficients,
%                           clean_poly after subs(x0.*x0,x0+1), adjoints on
%                           different output spaces (error), X.^2 = X.*X
%                           for X = [x0;x], and x.'*y = x'*y;
%   case_review_programs    piesos_quadvar with opts.L only, piesos_sosvar
%                           with exclude [1;0;0] and [1;1;1] and an ODE
%                           state, and degs with fields L and R with a
%                           dropped entry;
%   case_review_quadpoly    int of a quadPoly over 2 variables, subs of a
%                           name and a number at once, and of 2 names with
%                           1 value (error).
%
% INPUTS
% - only:   (optional) cellstr of names of the cases to run. Defaults to all.
%
% DJ (CC), 11/24/2026: Initial coding;

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'intop'),fullfile(here,'tensopvar'),fullfile(here,'polyopvar'));
addpath(genpath(fullfile(here,'..','..','sopvar')),genpath(fullfile(here,'..','..','misc')));

cases = {@case_eq_multirow, @case_liediff_vector, @case_intop_diff, @case_scalar_products, @case_constructor, @case_otimes_params, @case_mtimes_unchecked, @case_quadpoly_poly, @case_intop_plus, @case_exec_linear, @case_quadpoly_fixes, @case_dpvar_constant, @case_known_function, @case_unknown_options, @case_minor_fixes, @case_functional_product, @case_vector_scalar, ...
         @case_quadopvar_order, @case_apply_poly_checks, @case_review_operators, ...
         @case_review_polyopvar, @case_review_programs, @case_review_quadpoly};
if nargin>=1 && ~isempty(only)
    cases = cases(ismember(cellfun(@func2str,cases,'UniformOutput',false),only));
end
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
    end
end
if nfail>0
    error('Test_PIESOS_edge_cases: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_PIESOS_edge_cases: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_eq_multirow()
x = polyopvar('x',1,{'s'},[0,1]);
dpvar c1 c2
f = [(1-c1)*(x'*x); (2-c2)*(x'*x)];
prog = piesos_program(x,[c1;c2]);
prog = piesos_eq(prog,f);
evalc('prog_sol = piesos_solve(prog,struct(''simplify'',true));');
v = [double(piesos_getsol(prog_sol,c1)); double(piesos_getsol(prog_sol,c2))];
assert(norm(v-[1;2])<1e-6,'c = [%g;%g], expected [1;2].',v(1),v(2));
end

function case_liediff_vector()
pvar t s
x = pde_var('state',2,s,[0,1]);
PDE = [diff(x,t)==diff(x,s,2)+[1,2;0,1]*x-x.^3; subs(x,s,0)==0; subs(x,s,1)==0];
evalc('PIE = convert(PDE,''silent'');');
nm = PIE.space_pde.x.varname(PIE.pde2pie.x.idx(:,1));
xu = polyopvar(nm,PIE.space.x.dims,PIE.space.vars,PIE.space.dom,PIE.space.x.depmat);
dV = Liediff(xu'*xu,PIE);
u = PIE.T*PIE.x;
[is_eq,err] = polyopvar_equal(dV,2*(u'*PIE.f),1e-9);
assert(is_eq,'Liediff(<u,u>) differs from 2*<T*x,f> (error %g).',err);
end

function case_intop_diff()
x = polyopvar('x',1,{'s'},[0,1]);
V = x'*x;
K = V.C.ops{1};
dK = diff(K,{'s'},[0;1]);
for k=1:numel(K.params)
    P = K.params{k};    Q = dK.params{k};
    assert(isequal(size(Q),[2*size(P,1),size(P,2)]),'Wrong size of the parameters.');
    assert(full(max(max(abs(Q(1:size(P,1),:)-P))))==0,'Order 0 does not return the functional.');
    assert(nnz(Q(size(P,1)+1:end,:))==0,'Order 1 is not zero.');
end
end

function case_scalar_products()
pvar t s
x = pde_var(s,[0,1]);
PDE = [diff(x,t)==diff(x,s,2)+5*x-x*diff(x,s); subs(x,s,0)==0; subs(x,s,1)==0];
evalc('PIE = convert(PDE,''silent'');');
u = PIE.T*PIE.x;
[is_eq,err] = polyopvar_equal(2*u'*PIE.f,2*(u'*PIE.f),1e-12);
assert(is_eq,'2*u''*f differs from 2*(u''*f) (error %g).',err);
y = polyopvar('y',1,{'s'},[0,1]);
[is_eq,err] = polyopvar_equal(y'*(y*3),y'*(3*y),1e-12);
assert(is_eq,'y*3 differs from 3*y (error %g).',err);
% A decision variable times a functional is accepted, times a state not
dpvar c
assert(isa((y'*y)*c,'polyopvar') && isa(c*(y'*y),'polyopvar'),'(y''*y)*c or c*(y''*y) failed.');
for k=1:2
    msg = '';
    try
        if k==1,    y*c;    else,   c*y;    end %#ok<VUNUS>
    catch err
        msg = err.message;
    end
    assert(contains(msg,'is not supported'),'A decision variable times a state gave "%s".',msg);
end
% Also for a state without spatial variables; the weight form y'*c*y works
xo = polyopvar({'xo'},1,{'s'},[0,1],false);
msg = '';
try, c*xo; catch err, msg = err.message; end %#ok<VUNUS>
assert(contains(msg,'is not supported'),'c*xo gave "%s".',msg);
W = y'*c*y;
assert(all(cellfun(@(op) isa(op,'intop'),W.C.ops)),'y''*c*y is not a functional.');
end

function case_constructor()
threw = false;
try, polyopvar({'x','x'},[1;1],{'s'},[0,1]); catch, threw = true; end
assert(threw,'Duplicate state names were accepted.');
threw = false;
try, polyopvar('x',1,{'s'},[1,0]); catch, threw = true; end
assert(threw,'A domain [1,0] was accepted.');
end

function case_otimes_params()
pvar s1
Q = quadPoly.polynomial2quadPoly(1+s1,{'s1'},{});
S = struct('C',{{Q.C}},'Zs',{Q.Zs},'Zt',{Q.Zt},'ns',{Q.ns},'nt',{Q.nt},'dim',Q.dim);
A = otimes_params(S,2,[1,1]);
B = otimes_params(2,S,[1,1]);
assert(isequal(full(A.C{1}),full(B.C{1})),'otimes_params depends on the order of a double and a struct.');
end

function case_mtimes_unchecked()
rng(5011);
dims = [2,1; 1,2];
vars = struct('in',{{'sShared'}},'out',{{'sOut1','sShared'}});
depmat = struct('in',true(2,1),'out',logical([1,1;1,1]));
Bop = rand_tensopvar(dims,vars,[0,1],depmat,1,[2,2],1,0.5);
Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,ones(1,numel(Bop.vars.out)));
C1 = mtimes_tensopvar(Aop,Bop);
Aop.chkval = false;
C2 = mtimes_tensopvar(Aop,Bop);
assert(isequal(C1.params,C2.params),'The composition depends on chkval.');
end

function case_quadpoly_poly()
pvar s1
P = quadPoly.polynomial2quadPoly(1+s1,{'s1'},{});
R = quadPoly.quadPoly2polynomial(P*(s1+2));
err = double(max(abs(subs(R-(1+s1)*(s1+2),s1,linspace(0,1,7)))));
assert(err<1e-12,'quadPoly*polynomial is wrong (error %g).',err);
end

function case_intop_plus()
x = polyopvar('x',1,{'s'},[0,1]);
V = x'*x;
msg = '';
try, V.C.ops{1}+2; catch err, msg = err.message; end
assert(contains(msg,'can only be added'),'Unclear or no error for intop+2: "%s".',msg);
end

function case_quadpoly_fixes()
pvar s1 s2 s3 t1
q2p = @(Q) quadPoly.quadPoly2polynomial(Q);
p2q = @(p,ns,nt) quadPoly.polynomial2quadPoly(p,ns,nt);
maxdiff = @(a,b,v) max(arrayfun(@(p) max(max(abs(double(subs(a-b,v,p))))),linspace(0,1,5)));
% Concatenation with numeric inputs
P = p2q(1+s1^2,{'s1'},{});
F = 1+s1^2;
assert(maxdiff(q2p([P,2]),[F,2],s1)<1e-12,'[P,2] is wrong.');
assert(maxdiff(q2p([P;2]),[F;2],s1)<1e-12,'[P;2] is wrong.');
assert(maxdiff(q2p(blkdiag(P,2)),blkdiag(F,2),s1)<1e-12,'blkdiag(P,2) is wrong.');
% Integration over variables P does not depend on
Q = p2q(3*s2,{'s2'},{});
I1 = q2p(int(Q,{'s1'},0,2));                        % 2*3*s2
assert(maxdiff(I1,6*s2,s2)<1e-12,'Definite integral over an absent variable is wrong.');
I2 = q2p(int(Q,{'t1'}));                            % antiderivative t1*3*s2
assert(maxdiff(subs(I2,t1,0.5),1.5*s2,s2)<1e-12,'Antiderivative in an absent right variable is wrong.');
I3 = q2p(int(Q,{'s1','s2'},[0,0],[1,1]));           % int_0^1 int_0^1 3*s2 = 3/2
assert(abs(double(I3)-1.5)<1e-12,'Double integral with an absent variable is wrong.');
% Linear indexing with a vector
M = [1+s1, 2; s1^2, 3*s3];
PM = p2q(M,{'s1','s3'},{});
Mv = M(:);
for k = {[1,2],[1;3],[4,1,2]}
    kk = k{1};
    R = PM(kk);
    assert(isequal(size(R),size(M(kk))),'P(k) has size %s, expected %s.',mat2str(size(R)),mat2str(size(M(kk))));
    assert(maxdiff(subs(q2p(R),s3,0.3),subs(reshape(Mv(kk),size(M(kk))),s3,0.3),s1)<1e-12,'P(k) has the wrong entries.');
end
assert(isequal(size(PM(:)),[4,1]),'P(:) is not a column.');
Prow = PM(1,:);
assert(isequal(size(Prow([2;1])),[1,2]),'Linear indexing of a row quadPoly is not a row.');
end

function case_dpvar_constant()
x = polyopvar('x',1,{'s'},[0,1]);
V = x'*x;
dpvar c
% Constant terms are combined into one 'dpvar' coefficient
for F = {V - c, c - V, V - 1 - c, V - c - 1, 2*(V - c), V - c + c}
    nconst = sum(~any(F{1}.Z,2));
    assert(nconst<=1,'%d separate constant terms.',nconst);
end
s = evalc('display(V - c + 2)');
assert(contains(s,'2 - c') || contains(s,'- c + 2'),'Constant terms are not displayed: %s',s);
% max c s.t. int(u^4 - u^2) - c >= 0: the minimum of u^4-u^2 is -1/4
prog = piesos_program(x,c);
prog = piesos_setobj(prog,-c);
evalc('prog = piesos_ineq(prog,(x*x)''*(x*x) - V - c);');
evalc('sol = piesos_solve(prog,struct(''simplify'',true));');
cval = double(piesos_getsol(sol,c));
assert(abs(cval+0.25)<1e-5,'max c = %g, expected -0.25.',cval);
% Equality V - c = V - 2
prog = piesos_program(x,c);
prog = piesos_eq(prog,(V - c) - (V - 2));
evalc('sol = piesos_solve(prog,struct(''simplify'',true));');
cval = double(piesos_getsol(sol,c));
assert(abs(cval-2)<1e-6,'c = %g, expected 2.',cval);
% A decision variable is not added to a function of s
threw = false;
try, x + c; catch, threw = true; end %#ok<VUNUS>
assert(threw,'x + c (x a function of s) was accepted.');
end

function case_known_function()
pvar s
x = polyopvar('x',1,{'s'},[0,1]);
X0 = {1+s^2};
v1 = double(eval_polyopvar(x'*(x + s),X0));
v2 = double(eval_polyopvar(x'*(s*x),X0));
v3 = double(eval_polyopvar(x'*(x*s),X0));
e1 = double(int((1+s^2)*(1+s^2+s),s,0,1));
e2 = double(int(s*(1+s^2)^2,s,0,1));
assert(abs(v1-e1)<1e-10,'int x*(x+s) = %g, expected %g.',v1,e1);
assert(abs(v2-e2)<1e-10 && abs(v3-e2)<1e-10,'int x*s*x = %g, %g, expected %g.',v2,v3,e2);
threw = false;
pvar r
try, x + r; catch, threw = true; end %#ok<VUNUS>
assert(threw,'x(s) + r was accepted, although x does not depend on r.');
end

function case_unknown_options()
x = polyopvar('x',1,{'s'},[0,1]);
prog = piesos_program(x);
wstate = warning('on','PIESOS:unknownOption');
s1 = evalc('piesos_ineq(prog,x''*x,struct(''degz'',3));');
s2 = evalc('piesos_ineq(prog,x''*x,struct(''degs'',2,''psatz'',1));');
warning(wstate);
n1 = numel(strfind(s1,'unrecognized option'));
n2 = numel(strfind(s2,'unrecognized option'));
assert(n1==1 && n2==0,'%d warnings for a misspelled option (expected 1), %d for valid ones (expected 0).',n1,n2);
end

function case_minor_fixes()
pvar s
at = @(Y) double(subs(Y,s,0.5));
x = polyopvar('x',1,{'s'},[0,1]);
T = x.C.ops{1};                     % identity on x
xx = x*x;   Txx = xx.C.ops{1};      % pointwise product x.*x
v = {1+s^2};                        % x(0.5) = 1.25
assert(norm(at(apply_op([1;2].*T,v))-[1.25;2.5])<1e-12,'[1;2].*T is wrong.');
assert(norm(at(apply_op([1;2].*Txx,{1+s^2,1+s^2}))-[1.5625;3.125])<1e-12,'[1;2].*(x.*x) is wrong.');
assert(norm(at(apply_op([1 2;3 4].*repmat(T,2,2),{[1+s^2;1+s^2]}))-[3.75;8.75])<1e-12,'[1 2;3 4].*repmat(T,2,2) is wrong.');
assert(abs(at(apply_op(kron(2,T),v))-2.5)<1e-12,'kron(2,T) is wrong.');
V = x'*x;   K = V.C.ops{1};
assert(iscell(K(1,1).params),'K(1,1).params did not return the parameters.');
% piesos_eq(V - 1): the constant term gives the constraint 1 = 0
prog = piesos_program(x);
prog = piesos_eq(prog,V - 1);
evalc('sol = piesos_solve(prog,struct(''simplify'',true));');
assert(sol.solinfo.info.pinf || sol.solinfo.info.dinf,'piesos_eq(V-1) was not infeasible.');
% piesos_getsol of c*(1-V), with c = 3: 3*(1 - int (1+s^2)^2) = 3*(1-28/15)
dpvar c
prog = piesos_program(x,c);
prog = piesos_eq(prog,constant2polyopvar(c) - 3);
evalc('sol = piesos_solve(prog,struct(''simplify'',true));');
F = piesos_getsol(sol,c*(1-V));
assert(abs(double(apply_poly(F,v))-3*(1-28/15))<1e-6,'piesos_getsol(c*(1-V)) is wrong.');
% The largest c with 2*V + c*(1-V) = (2-c)*V + c >= 0 for all x, and c >= 0, is 2
prog = piesos_program(x,c);
prog = piesos_ineq(prog,c);
prog = piesos_ineq(prog,2*V + c*(1-V));
prog = piesos_setobj(prog,-c);
evalc('sol = piesos_solve(prog,struct(''simplify'',true));');
assert(abs(double(piesos_getsol(sol,c))-2)<1e-4,'The largest c with 2*V + c*(1-V) >= 0 is %g, expected 2.',double(piesos_getsol(sol,c)));
end

function case_functional_product()
pvar s
x = polyopvar('x',1,{'s'},[0,1]);
xo = polyopvar({'xo'},1,{'s'},[0,1],false);
V = x'*x; %#ok<NASGU>
cs = {'x*V','V*x','(1+V)*x','x*(1+V)','otimes(x,V)','otimes(V,x)'};
for k=1:numel(cs)
    msg = '';
    try, eval([cs{k} ';']); catch err, msg = err.message; end
    assert(contains(msg,'product of a functional'),'%s gave "%s".',cs{k},msg);
end
% A scalar ODE state times a functional is still supported
P = xo*(x'*x);
vals = cell(1,numel(P.varname));
for k=1:numel(P.varname)
    if strcmp(P.varname{k},'x'), vals{k} = 1+s^2; else, vals{k} = polynomial(2); end
end
v = double(apply_poly(P,vals));
assert(abs(v-2*28/15)<1e-10,'xo*(x''*x) = %g, expected %g.',v,2*28/15);
end

function case_vector_scalar()
pvar s
u = polyopvar('u',1,{'s'},[0,1]);
y = polyopvar('y',2,{'s'},[0,1]);
xo = polyopvar({'xo'},1,{'s'},[0,1],false);
V = u'*u; %#ok<NASGU>
vals_of = struct('u',1+s^2,'y',[1+s;3-s],'xo',polynomial(2));
% Values at s = 0.5: u = 1.25, y = [1.5;2.5], xo = 2, V = 28/15
U = 1.25;  Y = [1.5;2.5];  X = 2;  I = 28/15;
cs = {'[u;u]*u', [U^2;U^2];   'u*[u;u]', [U^2;U^2];   '[u;xo]*u', [U^2;X*U];
      'y*u', Y*U;             '[y;u]*u', [Y*U;U^2];   '[V;V]*V', [I^2;I^2]};
for k=1:size(cs,1)
    P = eval(cs{k,1});
    vals = cellfun(@(n) vals_of.(n),P.varname(:)','UniformOutput',false);
    v = double(subs(apply_poly(P,vals),s,0.5));
    assert(isequal(size(v),size(cs{k,2})) && max(abs(v-cs{k,2}))<1e-10, ...
        '%s = %s, expected %s.',cs{k,1},mat2str(v',6),mat2str(cs{k,2}',6));
end
% Vectors of equal length with different row-blocks: a clear error
msg = '';
try, y.*[u;u]; catch err, msg = err.message; end %#ok<VUNUS>
assert(contains(msg,'same row-blocks'),'y.*[u;u] gave "%s".',msg);
end

function case_quadopvar_order()
x1 = polyopvar('x1',polynomial({'s1';'s2'}),[0 1;0 1]);
prog = piesos_program(x1);
% Variables listed as (s2,s1), a monomial of 2 factors with tensor-product
% output (2 copies of each variable): intmat {'steps','full'} gives steps
% (2 regions) in s2 and the full domain in s1, as psatz [1 0] would act in s2
opts = struct('otype',2,'intmat',{{'steps','full'}});
[~,~,ZL,~,intmat] = piesos_quadopvar(prog,1,{'s2','s1'},[0 1;0 1],[2 2],0,opts);
assert(isequal(ZL.vars.out(:)',{'s1','s2'}),'Zop output variables are not (s1,s2).');
assert(size(intmat{1},1)==1 && size(intmat{2},1)==2,'intmat has %d regions in s1 and %d in s2, expected 1 and 2.', ...
    size(intmat{1},1),size(intmat{2},1));
% A region in a variable on which no state depends: clear error
msg = '';
try, piesos_quadopvar(prog,1,{'s1','s2'},[0 1;0 1],[2 0],0,struct('otype',2,'intmat',{{'full','steps'}})); catch err, msg = err.message; end
assert(contains(msg,'no state depends'),'intmat in an unused variable gave "%s".',msg);
end

function case_apply_poly_checks()
pvar s
x = polyopvar('x',1,{'s'},[0,1]);
y = polyopvar('y',1,{'s'},[0,1]);
P = otimes(x,y,2);                  % tensor-product output, depmat.out = 2
msg = '';
try, apply_poly(P,{1+s}); catch err, msg = err.message; end
assert(contains(msg,'Values must be specified'),'Missing value gave "%s".',msg);
dpvar c
msg = '';
try, apply_poly(x,{c*(1+s)}); catch err, msg = err.message; end
assert(contains(msg,'decision variables'),'A value with decision variables gave "%s".',msg);
% a value without decision variables still works
v = apply_poly(x,{polynomial(2)});
assert(abs(double(v)-2)<1e-12,'apply_poly(x,{2}) is wrong.');
end

function case_review_operators()
pvar s t
at = @(Y) double(subs(polynomial(Y),s,0.5));
% T*P for a 1x1 'sopvar' P on the right composes as T(P(x)), as T*tensopvar(P)
u = pde_var(s,[0,1]);
PDE = [diff(u,t)==diff(u,s,2); subs(u,s,0)==0; subs(u,s,1)==0];
evalc('PIE = convert(PDE,{''silent'',''nonlinear''});');
Tb = PIE.T.ops{1};
if isa(Tb,'sopvar'), Tb = tensopvar(Tb); end
P2 = mat2sopvar(s,struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0 1],'out',[0 1]));
v = {1+s};
ya = at(apply_op(Tb*P2,v));     yb = at(apply_op(Tb*tensopvar(P2),v));     yc = at(apply_op(Tb,{s*(1+s)}));
assert(abs(ya-yc)<1e-12 && abs(yb-yc)<1e-12,'T*P = %g, T*tensopvar(P) = %g, T(s*v) = %g.',ya,yb,yc);
% [1,2].*T and [1;2].*T for T = diag(1+s,2): both diag(1+s,4)
P = mat2sopvar([1+s 0;0 2],struct('in',{{'s'}},'out',{{'s'}}),struct('in',[0 1],'out',[0 1]));
T = tensopvar(P);
y1 = at(apply_op([1,2].*T,{[1;1]}));     y2 = at(apply_op([1;2].*T,{[1;1]}));
assert(norm(y1-[1.5;4])<1e-12 && norm(y2-[1.5;4])<1e-12,'[1,2].*T or [1;2].*T is wrong.');
% T{:,1}: the first factor of the (single) term
assert(isa(T{:,1},'sopvar'),'T{:,1} did not return the factor.');
% V.*[V;V*V] and [1;2]*(x.*x), for x = 1+s^2 (int x^2 = 28/15)
x = polyopvar('x',1,{'s'},[0,1]);
V = x'*x;   I = 28/15;   xv = {1+s^2};
W = V.*[V;V*V];
assert(norm(double(apply_poly(W,xv))-[I^2;I^3])<1e-10,'V.*[V;V*V] is wrong.');
assert(norm(at(apply_poly([1;2]*(x.*x),xv))-[1.25^2;2*1.25^2])<1e-12,'[1;2]*(x.*x) is wrong.');
% A functional times a 'tensopmat' with an empty block: an empty block
one = constant2polyopvar(1,{'s'},[0,1]);
L = one'*x;     K = L.C.ops{1};
B = [x.C, x.C];
B.ops{1,2} = [];
KB = K*B;
assert(isa(KB.ops{1,1},'intop') && isempty(KB.ops{1,2}),'K*[T,0] is wrong.');
end

function case_review_polyopvar()
pvar s t
x = polyopvar('x',1,{'s'},[0,1]);  y = polyopvar('y',1,{'s'},[0,1]);
Yv = 2-s;
% subs of a constant into int x*x*y: 4*int y (int (2-s) = 1.5)
W = subs(x'*(x.*y),{'x'},2);
assert(isequal(W.varname(:),{'y'}) && abs(double(apply_poly(W,{Yv}))-6)<1e-10,'subs(x''*(x.*y),{x},2) is not 4*int y.');
% apply_poly with 'dpvar' coefficients: c*int (1+s^2)^2 = 28/15*c
dpvar c
Y = apply_poly(c*(x'*x),{1+s^2});
D = Y - (28/15)*c;
assert(isa(Y,'dpvar') && full(max(abs(D.C(:))))<1e-10,'apply_poly(c*(x''*x)) is wrong.');
% clean_poly after subs(x0.*x0,x0+1): (x0+1)^2 = 9 at x0 = 2
x0 = polyopvar('x0');
Q = clean_poly(subs(x0.*x0,{'x0'},x0+1));
assert(abs(double(apply_poly(Q,{2}))-9)<1e-10,'clean_poly(subs(x0.*x0,x0+1)) is wrong.');
% Adjoints on different output spaces: a clear error
xs = polyopvar('x',1,{'s'},[0,2]); yt = polyopvar('y',1,{'t'},[0,2]); vs = polyopvar('v',1,{'s'},[0,2]);
msg = '';
try, (xs'+yt')*vs; catch err, msg = err.message; end %#ok<VUNUS>
assert(contains(msg,'same output space'),'(x''+y'')*v gave "%s".',msg);
% X.^2 = X.*X for X = [x0;x]
X = [x0;x];
D2 = X.^2 - X.*X;
vals = struct('x0',polynomial(3),'x',1+s);
r = apply_poly(D2,cellfun(@(n) vals.(n),D2.varname(:)','UniformOutput',false));
assert(max(abs(double(subs(polynomial(r),s,0.5))))<1e-10,'X.^2 differs from X.*X for X = [x0;x].');
% x.'*y = x'*y
assert(abs(double(apply_poly(x.'*y,{1+s,Yv}))-double(apply_poly(x'*y,{1+s,Yv})))<1e-12,'x.''*y differs from x''*y.');
end

function case_review_programs()
% opts.L alone: both bases otype 2, of the same size
x2 = polyopvar({'x1','x2'},[1;1],{'s'},[0,1]);
Z = dmonomials(x2,2,'full');
prog = piesos_program(x2);
[~,~,~,ZL,ZR] = piesos_quadvar(prog,Z,[],struct('L',struct('otype',2)));
oL = cellfun(@(K) K.otype(1),ZL.ops(~cellfun('isempty',ZL.ops)));
oR = cellfun(@(K) K.otype(1),ZR.ops(~cellfun('isempty',ZR.ops)));
assert(isequal(size(ZL),size(ZR)) && all(oL==2) && all(oR==2),'opts.L alone is not used for both bases.');
xo = polyopvar({'xo'},1,{'s'},[0,1],false);
% exclude [1;0;0] excludes the multiplier, the only term, of the ODE state
% (monopmat): V does not depend on it; all terms excluded is an error
x = polyopvar('x',1,{'s'},[0,1]);
prog = piesos_program([x;xo]);
[~,V] = piesos_sosvar(prog,[x;xo],struct('exclude',[1;0;0]));
assert(isequal(full(V.Z),[1,1]),'V depends on the ODE state with exclude [1;0;0].');
msg = '';
try, piesos_sosvar(prog,[x;xo],struct('exclude',[1;1;1])); catch err, msg = err.message; end
assert(contains(msg,'All terms'),'Excluding all terms gave "%s".',msg);
% degs as a struct with fields L and R (not per entry), with an entry
% dropped by opts.unique
x1 = polyopvar('x1',1,{'s'},[0,1]);
Zb = [otimes(x1,xo); otimes(xo,x1)];
prog = piesos_program([x1;xo]);
dg = struct('in',1,'out',1);
piesos_quadvar(prog,Zb,[],struct('degs',struct('L',dg,'R',dg)));
end

function case_review_quadpoly()
pvar s1 s2
q2p = @(Q) quadPoly.quadPoly2polynomial(Q);
Q = quadPoly.polynomial2quadPoly(polynomial(1),{'s1','s2'},{});
assert(abs(double(q2p(int(Q,{'s1','s2'},[1,1],[2,2])))-1)<1e-12,'The integral of 1 over [1,2]^2 is not 1.');
% Iterated: int_0^1 int_0^{s2} 1 ds1 ds2 = 1/2
assert(abs(double(q2p(int(Q,{'s1','s2'},{0,0},{'s2',1})))-0.5)<1e-12,'The iterated integral is not 1/2.');
P = quadPoly.polynomial2quadPoly(s1+2*s2^2,{'s1','s2'},{});
R = q2p(subs(P,{'s1','s2'},{'s2',0}));
assert(abs(double(subs(R,s2,0.7))-0.7)<1e-12,'subs(P,{s1,s2},{s2,0}) is not P(s2,0).');
msg = '';
try, subs(P,{'s1','s2'},{1}); catch err, msg = err.message; end
assert(contains(msg,'number of vars and vals'),'subs with 2 names and 1 value gave "%s".',msg);
end

function case_exec_linear()
pvar t s
st = PIESOS_settings('light');     st.verbose = false;
r_vals = [5,10.5];      is_stable = [true,false];
for k=1:2
    x = pde_var(s,[0,1]);
    PDE = [diff(x,t)==diff(x,s,2)+r_vals(k)*x; subs(x,s,0)==0; subs(x,s,1)==0];
    evalc('PIE = convert(PDE,''silent'');');
    assert(isa(PIE,'pie_struct'),'The linear PDE did not give a pie_struct.');
    evalc('[~,info] = PIESOS_stability_global(PIE,0,''pie'',st);');
    assert(info.feasible==is_stable(k),'r = %g: feasible = %d, expected %d.',r_vals(k),info.feasible,is_stable(k));
end
end
