% Test_polyopvar_subs_const_partial.m -- substitution of only some of the
% states in a monomial by constants ('double' or 'polynomial' values,
% subs_const in @polyopvar/subs.m). The factors replaced remain in the
% operator as multipliers from R (zeros in Z, see polyopvar) with the
% constant as value, so that e.g. x*y with x=c gives Z=[0,y]. Parts:
%   (A) two factors, pointwise and Kronecker output, by name and by
%       position;
%   (B) three factors, one or two of them replaced;
%   (C) a state repeated in the monomial, a o b o a;
%   (D) states depending on different spatial variables (and on none);
%   (E) vector-valued states;
%   (F) several monomials, with unaffected, partial and full ones, and
%       (F2) combination of consecutive multipliers from R after partial
%       substitution (pointwise output);
%   (G) functionals (intop), and further operations on the result;
%   (I) subs_var, by name, with a polyopvar for states depending on some
%       of the variables only;
%   (H) errors.
%
% Oracle: the substituted polynomial is evaluated at the remaining states
% (eval_polyopvar) and compared with the original polynomial evaluated
% at the constants and the same remaining states.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Evaluation uses eval_polyopvar (Testfolder/polyopvar), not
%   apply_poly: apply_poly evaluates by subs, which this test checks.
% - Expected Z for a pointwise ('pw') output has consecutive multipliers
%   from R combined (compress_zeros); for a Kronecker output they are
%   left as is.

pvar s1 s2
svar_all = {s1,s2};
tol = 1e-8;
ncases = 0;
allnames = {'a','b','c'};

%% (A)-(C) states depending on all variables, scalar
for N = 1:2
    svs = svar_all(1:N);
    dep = true(3,N);        sz = [1;1;1];
    for kind = {'pw','kron'}
        for form = {'double','polynomial'}
            rng(100*N+10*numel(kind{1})+numel(form{1}));
            lbl = sprintf('N=%d,%s,%s',N,kind{1},form{1});
            xv = struct('a',rand_const(svs,1,'polynomial'),'b',rand_const(svs,1,'polynomial'),'c',rand_const(svs,1,'polynomial'));
            ca = rand_const(svs,1,form{1});     cb = rand_const(svs,1,form{1});     cc = rand_const(svs,1,form{1});
            % (A) a o b
            f = make_f(N,[1,2],dep,sz,kind{1});
            ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2],tol,lbl+" (A) a");
            ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0],tol,lbl+" (A) b");
            % (B) a o b o c
            f = make_f(N,[1,2,3],dep,sz,kind{1});
            ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,3],tol,lbl+" (B) a");
            ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,3],tol,lbl+" (B) b");
            ncases = ncases+check_partial(f,{'c'},{cc},xv,sz,[1,2,0],tol,lbl+" (B) c");
            ncases = ncases+check_partial(f,{'a';'c'},{ca;cc},xv,sz,[0,2,0],tol,lbl+" (B) a,c");
            ncases = ncases+check_partial(f,{'c';'a'},{cc;ca},xv,sz,[0,2,0],tol,lbl+" (B) c,a");
            % (C) a o b o a: the state in two factors
            f = make_f(N,[1,2,1],dep,sz,kind{1});
            ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,0],tol,lbl+" (C) a");
            ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,1],tol,lbl+" (C) b");
            % one of the two factors a, by position
            % (evaluation by state gives all factors of a state the same value: apply the operator to the factors)
            Yref = apply_op(f.C,{{xv.a,xv.b,ca}});
            fnew = subs(f,3,ca);
            assert_close(eval_new(fnew,xv,sz,allnames),Yref,tol,lbl+" (C) factor 3");
            if ~isequal(pad(full(fnew.Z),3),[1,2,0]),   error("Test_polyopvar_subs_const_partial failed (%s): Z = [%s].",lbl,num2str(full(fnew.Z))),   end
            ncases = ncases+1;
        end
    end
end

%% (D) states depending on different variables: a(s1,s2), b(s1), c(no variables)
dep = [true,true; true,false; false,false];     sz = [1;1;1];
N = 2;      svs = svar_all;
for kind = {'pw','kron'}
    for form = {'double','polynomial'}
        rng(500+10*numel(kind{1})+numel(form{1}));
        lbl = sprintf('(D) %s,%s',kind{1},form{1});
        xv = struct('a',rand_const(svs,1,'polynomial'),'b',rand_const({s1},1,'polynomial'),'c',rand_const({},1,'polynomial'));
        ca = rand_const(svs,1,form{1});     cb = rand_const({s1},1,form{1});    cc = rand_const({},1,form{1});
        f = make_f(N,[1,2,3],dep,sz,kind{1});
        ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,3],tol,lbl+" a");
        ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,3],tol,lbl+" b");
        ncases = ncases+check_partial(f,{'c'},{cc},xv,sz,[1,2,0],tol,lbl+" c");
        ncases = ncases+check_partial(f,{'b';'c'},{cb;cc},xv,sz,[1,0,0],tol,lbl+" b,c");
        ncases = ncases+check_partial(f,{'a';'b'},{ca;cb},xv,sz,[0,0,3],tol,lbl+" a,b");
    end
end

%% (E) vector-valued states, sizes [2;3;1]
sz = [2;3;1];
for N = 1:2
    svs = svar_all(1:N);    dep = true(3,N);
    for form = {'double','polynomial'}
        rng(700+N+numel(form{1}));
        lbl = sprintf('(E) N=%d,%s',N,form{1});
        xv = struct('a',rand_const(svs,2,'polynomial'),'b',rand_const(svs,3,'polynomial'),'c',rand_const(svs,1,'polynomial'));
        ca = rand_const(svs,2,form{1});     cb = rand_const(svs,3,form{1});     cc = rand_const(svs,1,form{1});
        f = make_f(N,[1,2,3],dep,sz,'pw');
        ncases = ncases+check_partial(f,{'a'},{ca},xv,sz,[0,2,3],tol,lbl+" a");
        ncases = ncases+check_partial(f,{'b'},{cb},xv,sz,[1,0,3],tol,lbl+" b");
        ncases = ncases+check_partial(f,{'c';'a'},{cc;ca},xv,sz,[0,2,0],tol,lbl+" c,a");
    end
end

%% (F) several monomials: a o b (partial), a (full), b (unaffected), plus a o b o c
N = 1;      svs = svar_all(1);      dep = true(3,N);    sz = [1;1;1];
rng(900);
xv = struct('a',rand_const(svs,1,'polynomial'),'b',rand_const(svs,1,'polynomial'),'c',rand_const(svs,1,'polynomial'));
fm = make_f(N,[1,2],dep,sz,'pw') + make_f(N,1,dep,sz,'pw') + make_f(N,2,dep,sz,'pw') + make_f(N,[1,2,3],dep,sz,'pw');
for form = {'double','polynomial'}
    ca = rand_const(svs,1,form{1});     cb = rand_const(svs,1,form{1});
    lbl = "(F) "+string(form{1});
    ncases = ncases+check_partial(fm,{'a'},{ca},xv,sz,[],tol,lbl+" a");
    ncases = ncases+check_partial(fm,{'b'},{cb},xv,sz,[],tol,lbl+" b");
    ncases = ncases+check_partial(fm,{'a';'b'},{ca;cb},xv,sz,[],tol,lbl+" a,b");
    % a state that appears in no monomial leaves the polynomial unchanged
    cz = rand_const(svs,1,form{1});
    fnew = subs(fm,{'c'},cz);
    xs = xv;    xs.c = cz;
    assert_close(eval_new(fnew,rmfield(xv,'c'),sz,allnames),eval_polyopvar(fm,{xs.a;xs.b;xs.c}),tol,lbl+" c");
    ncases = ncases+1;
end

%% (F2) After partial substitution, consecutive multipliers from R are combined (pointwise output)
% and the result is cleaned: monomials that coincide are added, states that no longer appear removed.
% a o b o c and b o a o c, with a,b replaced: both Z=[0,0,c] --> [0,c] and then one monomial
N = 1;      svs = svar_all(1);      dep = true(3,N);    sz = [1;1;1];
for kind = {'pw','kron'}
    rng(950+numel(kind{1}));
    lbl = "(F2) "+string(kind{1});
    xv = struct('a',rand_const(svs,1,'polynomial'),'b',rand_const(svs,1,'polynomial'),'c',rand_const(svs,1,'polynomial'));
    ca = rand_const(svs,1,'double');        cb = rand_const(svs,1,'polynomial');
    fpair = make_f(N,[1,2,3],dep,sz,kind{1}) + make_f(N,[2,1,3],dep,sz,kind{1});
    xs = xv;    xs.a = ca;      xs.b = cb;
    Yref = eval_polyopvar(fpair,{xs.a;xs.b;xs.c});
    fnew = subs(fpair,{'a';'b'},[ca;cb]);
    assert_close(eval_new(fnew,rmfield(xv,{'a','b'}),sz,allnames),Yref,tol,lbl);
    Zexp = [0,0,1];     if strcmp(kind{1},'pw'),    Zexp = [0,1];   end
    n = max(size(fnew.Z,2),numel(Zexp));
    if size(fnew.Z,1)~=1 || ~isequal(pad(full(fnew.Z),n),pad(Zexp,n)) || any(ismember({'a','b'},fnew.varname))
        error("Test_polyopvar_subs_const_partial failed (%s): Z = [%s] in states {%s}, expected the single monomial [%s] in state c.",lbl,num2str(full(fnew.Z)),strjoin(fnew.varname,','),num2str(Zexp))
    end
    ncases = ncases+2;
end

%% (G) functionals, and further operations on the result
polyopvar a(s1) b(s1)
V = innerprod(a,b);         % int a*b ds, Z=[1,2]
rng(1000);
xa = rand_const({s1},1,'polynomial');       xb = rand_const({s1},1,'polynomial');
vals = struct('a',xa,'b',xb);
for form = {'double','polynomial'}
    cc = rand_const({s1},1,form{1});
    tests = {   subs(V,{'b'},cc),       int(xa*cc,s1,0,1)
                subs(V,2,cc),           int(xa*cc,s1,0,1)
                subs(V,{'a'},cc),       int(cc*xb,s1,0,1)       };
    for k=1:size(tests,1)
        assert_close(eval_V(tests{k,1},vals),tests{k,2},tol,sprintf('(G) functional %d, %s',k,form{1}));
        ncases = ncases+1;
    end
end
% further operations: product and sum with other polynomials of the same output space
N = 1;      svs = svar_all(1);      dep = true(3,N);    sz = [1;1;1];
rng(1100);
xv = struct('a',rand_const(svs,1,'polynomial'),'b',rand_const(svs,1,'polynomial'),'c',rand_const(svs,1,'polynomial'));
ca = rand_const(svs,1,'double');
f = make_f(N,[1,2],dep,sz,'pw');        fz = make_f(N,3,dep,sz,'pw');
g = subs(f,{'a'},ca);
xs = xv;    xs.a = ca;
Yf = eval_polyopvar(f,{xs.a;xs.b;xs.c});       Yz = eval_polyopvar(fz,{xs.a;xs.b;xs.c});
assert_close(eval_new(g.*fz,rmfield(xv,'a'),sz,allnames),Yf.*Yz,tol,'(G) times');
assert_close(eval_new(fz.*g,rmfield(xv,'a'),sz,allnames),Yz.*Yf,tol,'(G) times, reversed');
assert_close(eval_new(g+fz,rmfield(xv,'a'),sz,allnames),Yf+Yz,tol,'(G) plus');
ncases = ncases+3;

%% (I) subs_var, by name, with a polyopvar for states depending on only some of the variables:
% a(s1,s2), b(s1), c(no variables). The output variables of the substitute are those of the
% state, a subset of those of f, so that the function spaces must be compared by variable name
N = 2;      svs = svar_all;     dep = [true,true; true,false; false,false];     sz = [1;1;1];
polyopvar w(s1) v(s1)
for kind = {'pw','kron'}
    rng(1200+numel(kind{1}));
    lbl = "(I) "+string(kind{1});
    xv = struct('a',rand_const(svs,1,'polynomial'),'b',rand_const({s1},1,'polynomial'),'c',rand_const({},1,'polynomial'));
    ca = rand_const(svs,1,'polynomial');    cb = rand_const({s1},1,'polynomial');   cc = rand_const({},1,'polynomial');
    ga = known_fun(ca,{'s1','s2'});         gb = known_fun(cb,{'s1'});              gc = known_fun(cc,cell(1,0));
    f = make_f(N,[1,2,3],dep,sz,kind{1});
    % known functions
    ncases = ncases+check_known(f,{'a'},{ga},{ca},xv,sz,tol,lbl+" a");
    ncases = ncases+check_known(f,{'b'},{gb},{cb},xv,sz,tol,lbl+" b");
    ncases = ncases+check_known(f,{'c'},{gc},{cc},xv,sz,tol,lbl+" c");
    ncases = ncases+check_known(f,{'b';'c'},{gb;gc},{cb;cc},xv,sz,tol,lbl+" b,c");
    ncases = ncases+check_known(f,{'a';'b'},{ga;gb},{ca;cb},xv,sz,tol,lbl+" a,b");
    % functions of new states w,v on the same variables as b
    xw = rand_const({s1},1,'polynomial');   xvv = rand_const({s1},1,'polynomial');
    for kv=1:2
        if kv==1,   val = 2*w;          Yb = 2*xw;
        else,       val = w+2*v;        Yb = xw+2*xvv;      end
        fnew = subs(f,{'b'},val);
        Yref = eval_polyopvar(f,{xv.a;Yb;xv.c});
        xs = struct('a',xv.a,'c',xv.c,'w',xw,'v',xvv);
        assert_close(eval_new(fnew,xs,sz,allnames),Yref,tol,lbl+sprintf(" b=val%d",kv));
        ncases = ncases+1;
    end
end
% a substitute on other variables than the state: b(s1) by a function of s2, a(s1,s2) by a function of s1
% only, or by a function of a variable that f does not have
f = make_f(N,[1,2,3],dep,sz,'pw');
expect_error(@() subs(f,{'b'},known_fun(1+s2,{'s2'})),'Output space');
expect_error(@() subs(f,{'a'},known_fun(1+s1,{'s1'})),'Output space');
expect_error(@() subs(f,{'b'},known_fun(1+s1,{'s3'})),'Output space');
ncases = ncases+3;

%% (H) errors
N = 2;      svs = svar_all;     dep = [true,true; true,false; false,false];     sz = [1;1;1];
f = make_f(N,[1,2,3],dep,sz,'pw');
% polynomial in a variable that the state does not depend on
expect_error(@() subs(f,{'b'},s2),'spatial variables');
expect_error(@() subs(f,{'c'},s1),'spatial variables');
expect_error(@() subs(f,2,s2),'spatial variables');
% wrong number of rows, unknown state, factor that is no state, factor beyond the last one
expect_error(@() subs(f,{'a'},[1;2]),'Number of rows');
expect_error(@() subs(f,{'d'},1),'states');
expect_error(@() subs(f,5,1),'degree');
% all states replaced in a function with a Kronecker output is still not supported
fk = make_f(N,[1,2],dep,sz,'kron');
expect_error(@() subs(fk,{'a';'b'},[1;2]),'tensor product');
ncases = ncases+7;

fprintf('Test_polyopvar_subs_const_partial passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function f = make_f(N,Zrow,dep,sz,kind)
% Function f = Kop*(x_{Z1} o ... o x_{Zd}) in the states a,b,c of sizes sz, depending on the
% variables of dep (3 x N, row i the variables of state i), with Kop random, with pointwise
% ('pw', output size 2) or Kronecker ('kron') output
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
d = numel(Zrow);
din = reshape(sz(Zrow),[],1);
if strcmp(kind,'pw')
    dims = [2*ones(d,1),din];       ot = [1,2];     mout = 2;
else
    dims = [(1+(1:d)'),din];        ot = [2,2];     mout = prod(dims(:,1));
end
dp = struct('in',dep(Zrow,:),'out',true(d,N));
Kop = rand_tensopvar(dims,vars,dom,dp,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = Zrow;                 f.varname = {'a';'b';'c'};
f.dims.in = sz(:);          f.dims.out = mout;
f.vars.in = svars;          f.dom.in = dom.in;          f.depmat.in = dep;
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
f = initialize(f);
end

function v = rand_const(svs,m,form)
% Random m x 1 constant function: 'double', or 'polynomial' in the variables svs (cell of pvar)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = polynomial(randn(m,1));
    for i=1:numel(svs)
        v = v + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
    end
end
end

function z = pad(z,n)
z(end+1:n) = 0;
end

function z = compress_zeros(z)
% Replace each run of consecutive zeros by a single zero
keep = true(size(z));
keep(2:end) = ~(z(2:end)==0 & z(1:end-1)==0);
z = z(keep);
end

function Y = eval_new(f,xvals,sz,allnames)
% Evaluate f at the values xvals.(name) of its states; states without a value get zeros
xv = cell(numel(f.varname),1);
for k=1:numel(xv)
    nm = f.varname{k};
    if isfield(xvals,nm)
        xv{k} = xvals.(nm);
    else
        xv{k} = polynomial(zeros(sz(strcmp(allnames,nm)),1));
    end
end
Y = eval_polyopvar(f,xv);
end

function n = check_partial(f,subnames,cvals,xvals,sz,Zexp,tol,label)
% Substitute the states subnames of f (by name, and by position if the monomial is a single one
% in which they each appear once) by the constants cvals, and compare with f evaluated at the
% constants. Zexp: expected Z of the result, if f has a single monomial (otherwise []).
allnames = {'a','b','c'};
cmat = vertcat(cvals{:});
xs = xvals;
for k=1:numel(subnames),    xs.(subnames{k}) = cvals{k};    end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
xrem = rmfield(xvals,subnames);
fnew = subs(f,subnames,cmat);
assert_close(eval_new(fnew,xrem,sz,allnames),Yref,tol,label+" (by name)");
if ~isempty(Zexp)
    % Expected Z: consecutive multipliers from R combined, for an operator with pointwise
    % output (not supported otherwise, so left), in the states of fnew, from which the states
    % that no longer appear have been removed (by name: Zexp is in terms of the states a,b,c)
    op = f.C.ops{1,1};
    if isa(op,'tensopvar') && op.otype(1)==1,   Zexp = compress_zeros(Zexp);    end
    for q=find(Zexp)
        idx = find(strcmp(fnew.varname,allnames{Zexp(q)}));
        if isempty(idx),    error("Test_polyopvar_subs_const_partial failed (%s): state %s was removed.",label,allnames{Zexp(q)});    end
        Zexp(q) = idx;
    end
    n = max(size(fnew.Z,2),numel(Zexp));
    if size(fnew.Z,1)~=1 || ~isequal(pad(full(fnew.Z),n),pad(Zexp,n))
        error("Test_polyopvar_subs_const_partial failed (%s): Z = [%s], expected [%s].",label,num2str(full(fnew.Z)),num2str(Zexp))
    end
end
n = 1;
Zr = full(f.Z);
idx = cellfun(@(s) find(strcmp(allnames,s)),subnames);
if size(Zr,1)==1 && all(arrayfun(@(i) sum(Zr==i)==1,idx))
    pos = arrayfun(@(i) find(Zr==i),idx);
    fnew = subs(f,pos(:)',cmat);
    assert_close(eval_new(fnew,xrem,sz,allnames),Yref,tol,label+" (by position)");
    n = 2;
end
end

function g = known_fun(m,vars)
% Degree-0 polyopvar representing the known function m (polynomial) of the variables vars
% (cellstr, in [0,1]): a multiplier from R to L2 on the dummy input 1
M = mat2tensopvar(m,struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
g = polyopvar();
g.varname = cell(0,1);      g.Z = sparse(1,0);      g.C = tensopmat(M);
g = initialize(g);
end

function n = check_known(f,subnames,gs,cvals,xvals,sz,tol,label)
% Substitute the states subnames of f by name by the known functions gs (polyopvar, stacked), and
% compare with f evaluated at the constant functions cvals that these represent
allnames = {'a','b','c'};
xs = xvals;
for k=1:numel(subnames),    xs.(subnames{k}) = cvals{k};    end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
fnew = subs(f,subnames,vertcat(gs{:}));
assert_close(eval_new(fnew,rmfield(xvals,subnames),sz,allnames),Yref,tol,label);
n = 1;
end

function Y = eval_V(V,vals)
% Evaluate the functional V (intop coefficients) at the states vals.(name); the dummy
% variables of s1 are numbered over the factors, as in apply_op of 'intop'
Y = 0;
Zf = full(V.Z);
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    v = 1;
    for j=1:numel(idcs)
        xj = vals.(V.varname{idcs(j)});
        if numel(idcs)>1,   xj.varname(strcmp(xj.varname,'s1')) = {['s1_',num2str(j)]};     end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function assert_close(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_subs_const_partial (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_subs_const_partial failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_polyopvar_subs_const_partial failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_polyopvar_subs_const_partial failed: expected an error containing '%s', but none was thrown.",str)
end
