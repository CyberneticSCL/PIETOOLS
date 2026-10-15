% Test_polyopvar_subs_class.m -- tests of @polyopvar/subs.m (all routes:
% subs_var, subs_factor, subs_const), complementing the
% Test_polyopvar_subs_* tests with a randomized sweep and the error and
% edge cases.
%
% Expected values are computed independently of subs: the original
% polynomial is evaluated (eval_polyopvar, or apply_op of the coefficient
% operators on the factors) at the constant or at the polynomial val(y)
% substituted, and compared with the result of the substitution evaluated
% at the remaining states.
% Parts:
%   (1) dispatch and input validation;
%   (2) renaming of states (subs_var, var and val cellstr);
%   (3) errors and edge cases of subs_var, subs_factor, subs_const;
%   (4) result is a valid polynomial with the output space of f
%       unchanged;
%   (5) sweep over 1-3 spatial variables, pointwise/Kronecker output,
%       vector-valued states, repeated factors: the first state and all
%       states, in random order, by name by a constant ('double',
%       'polynomial', 'dpvar') and by a function of a new state (sorted
%       between the states, so that the states are renumbered), the
%       first and last factor by position, and (5b) polynomials with
%       several monomials and several rows;
%   (6) functionals (intop coefficients).
% Cost: about 1.5 minutes.
%
% DJ, 10/15/2026: Initial coding;

tol = 1e-8;
ncases = 0;
pvar s1 s2 s3
svar_all = {s1,s2,s3};
allnames = {'a','b','c'};
sz1 = [1;1;1];

%% (1) dispatch and input validation
rng(1);
f = make_f(1,[1,2,3],sz1,'pw');
xv = struct('a',rand_pol(svar_all(1),1),'b',rand_pol(svar_all(1),1),'c',rand_pol(svar_all(1),1));
yy = polyopvar('y',{'s1'},[0,1]);
expect_error(@() subs(5,{'a'},yy),"must be of type 'polyopvar'");        % dispatches to polyopvar/subs through val
expect_error(@() subs(f,'a',1),"must be expressed as 'cellstr'");
expect_error(@() subs(f,'a',yy),"must be expressed as 'polyopvar' or numeric");
expect_error(@() subs(f,{'a'},'y'),'');
ncases = ncases+4;
% f that has not been validated is initialized first
f2 = f;     f2.chkval = false;
c = randn;
assert_close(eval_new(subs(f2,{'a'},c),xv,sz1,allnames),eval_polyopvar(f,{polynomial(c);xv.b;xv.c}),tol,'(1) chkval=false');
% ... also if its fields are as set by the user, not yet validated by initialize
rng(1);      g = make_f(1,[1,2,3],sz1,'pw',true);       % same operator as f
assert_close(eval_new(subs(g,{'a'},c),xv,sz1,allnames),eval_polyopvar(f,{polynomial(c);xv.b;xv.c}),tol,'(1) not initialized');
ncases = ncases+2;

%% (2) renaming of states
rng(2);
f = make_f(1,[1,2,3],sz1,'pw') + make_f(1,[1,3,2],sz1,'pw');
Yref = eval_polyopvar(f,{xv.a;xv.b;xv.c});
% a --> z moves the state after c; names given as column or row, states not in f are ignored
fnew = subs(f,{'a';'q'},{'z';'w'});
assert(isequal(fnew.varname,{'b';'c';'z'}),'(2) rename: states not sorted');
assert_close(eval_polyopvar(fnew,{xv.b;xv.c;xv.a}),Yref,tol,'(2) a --> z');
fnew = subs(f,{'a','b'},{'b','a'});
assert_close(eval_polyopvar(fnew,{xv.b;xv.a;xv.c}),Yref,tol,'(2) swap a,b');
% names listed in another order than the states
fnew = subs(f,{'c','a'},{'x','y'});
assert(isequal(fnew.varname,{'b';'x';'y'}),'(2) rename, var not in the order of the states');
assert_close(eval_polyopvar(fnew,{xv.b;xv.c;xv.a}),Yref,tol,'(2) c,a --> x,y');
fnew = subs(f,{'q'},{'w'});
assert(isequal(fnew.varname,f.varname) && isequal(fnew.Z,f.Z),'(2) rename of a state that is not in f');
% c --> b merges the states, and the monomials a*b*b that now coincide
fnew = subs(f,{'c'},{'b'});
assert(isequal(fnew.varname,{'a';'b'}) && size(fnew.Z,1)==1,'(2) c --> b: states or monomials not merged');
assert_close(eval_polyopvar(fnew,{xv.a;xv.b}),eval_polyopvar(f,{xv.a;xv.b;xv.b}),tol,'(2) c --> b');
% states of different size cannot be merged
fv = make_f(1,[1,2],[2;1;1],'pw');
expect_error(@() subs(fv,{'a'},{'b'}),'same function space');
ncases = ncases+8;

%% (3) errors and edge cases
rng(3);
svars = {'s1'};     dom = struct('in',[0,1],'out',[0,1]);
f = make_f(1,[1,2,3],sz1,'pw');
val1 = make_val(svars,dom,'y',1,1,1);
val2 = [val1;val1];
% subs_var
expect_error(@() subs(f,yy+yy,val1),'individial states');
expect_error(@() subs(f,[yy;yy],val2),'unique');
expect_error(@() subs(f,{'a','b'},val1),'Number of rows');
expect_error(@() subs(f,{'a'},make_val(svars,dom,'y',2,1,1)),'Output space');
expect_error(@() subs(f,{'a'},make_val({'s2'},struct('in',[0,1],'out',[0,1]),'y',1,1,1)),'Output space');
% names that are not states of f are ignored, together with their rows of val: nothing is substituted
fnew = subs(f,{'q'},val1);
assert(isequal(fnew.Z,f.Z) && isequal(fnew.varname,f.varname),'(3) unknown state should leave f unchanged');
vq = make_val(svars,dom,'y',1,1,1);     va = make_val(svars,dom,'y',1,1,2);
yv = struct('b',xv.b,'c',xv.c,'y',rand_pol(svar_all(1),1));
Yref = eval_new(subs(f,{'a'},va),yv,[sz1;1],[allnames,{'y'}]);
assert_close(eval_new(subs(f,{'q','a'},[vq;va]),yv,[sz1;1],[allnames,{'y'}]),Yref,tol,'(3) unknown and known state');
expect_error(@() subs(f,{'a','q','r'},val2),'Number of rows');ncases = ncases+6;
% subs_factor
expect_error(@() subs(f,1.5,val1),'unique positive integers');
expect_error(@() subs(f,0,val1),'unique positive integers');
expect_error(@() subs(f,-1,val1),'unique positive integers');
expect_error(@() subs(f,[1,1],val2),'unique positive integers');
expect_error(@() subs(f,4,val1),'cannot exceed');
expect_error(@() subs(f,1,'y'),'non-''polyopvar''');
expect_error(@() subs(f,[1,2],val1),'Number of rows');
expect_error(@() subs(f,1,make_val(svars,dom,'y',2,1,1)),'Input function space');       % size of the state
expect_error(@() subs(f,1,make_val({'s2'},dom,'y',1,1,1)),'Input function space');            % variable of f
expect_error(@() subs(make_f(2,[1,2],sz1,'pw'),1,make_val({'s1'},dom,'y',1,1,1)),'Input function space');    % fewer variables than the state
expect_error(@() subs(f+make_f(1,[1,2],sz1,'pw'),1,val1),'multiple monomials');
fk = subs(f,{'a'},randn);               % Z=[0,2,3]: the factor 1 is a multiplier from R
expect_error(@() subs(fk,1,val1),'not multipliers');
ncases = ncases+9;
% subs_const
dpvar dd
expect_error(@() subs(f,{'a'},dd),'decision variables');
expect_error(@() subs(f,{'a','b'},[1,2]),'column vector');
expect_error(@() subs(f,yy+yy,1),'individial states');
expect_error(@() subs(f,{'q'},1),'must be states');
expect_error(@() subs(f,{'a','a'},[1;2]),'unique');
expect_error(@() subs(f,{'a'},[1;2]),'Number of rows of constant');
expect_error(@() subs(f,{'a'},s2),'spatial variables');
expect_error(@() subs(f,{'a','b','c'},[1;2;s2]),'spatial variables');
expect_error(@() subs(f+make_f(1,[1,2],sz1,'pw'),1,1),'single monomial');
expect_error(@() subs(f,[1,4],[1;2]),'cannot exceed|not exceeding');
expect_error(@() subs(f,[1,1],[1;2]),'unique');
expect_error(@() subs(fk,1,1),'not multipliers');
expect_error(@() subs(f,[1,2],1),'Number of rows of constant');
expect_error(@() subs(make_f(1,[1,2],sz1,'kron'),{'a','b'},[1;2]),'tensor product');
% a state that is not in any monomial leaves f unchanged
fnew = subs(make_f(1,[1,2],sz1,'pw'),{'c'},1);
assert(isequal(size(fnew.Z),[1,2]),'(3) unused state should leave f unchanged');
% dpvar constants without decision variables are accepted
cc = randn;
assert_close(eval_new(subs(f,{'a'},dpvar(cc)),xv,sz1,allnames),eval_polyopvar(f,{polynomial(cc);xv.b;xv.c}),tol,'(3) dpvar constant');
ncases = ncases+16;

%% (4) the result is a valid polynomial, and keeps the output space of f
rng(4);
for N = 1:2
    [svars,dom] = spatial_space(N);
    for kind = {'pw','kron'}
        f = make_f(N,[1,2,3],sz1,kind{1});
        val = make_val(svars,dom,'y',1,1,2);
        cases = {subs(f,{'a'},val), subs(f,2,val), subs(f,{'a'},randn), subs(f,2,randn), subs(f,{'a'},{'z'})};
        for k=1:numel(cases)
            fnew = cases{k};
            g = fnew;   g.chkval = false;   g = initialize(g);
            assert(isequal(fnew.vars.out,f.vars.out) && isequal(fnew.dom.out,f.dom.out) && isequal(fnew.dims.out,f.dims.out) ...
                && isequal(fnew.depmat.out,f.depmat.out) && fnew.chkval && isequal(g.Z,fnew.Z),sprintf('(4) N=%d,%s,%d',N,kind{1},k));
            ncases = ncases+1;
        end
    end
end

%% (5) randomized sweep: functions
Zlist = {[1,2],[2,3],[1,2,3],[1,2,1],[3,1,1]};
for N = 1:3
    [svars,dom] = spatial_space(N);     svs = svar_all(1:N);
    for kind = {'pw','kron'}
        for iz = 1:numel(Zlist)
            Z = Zlist{iz};
            % cost: skip 3 factors in 2D with Kronecker output, and 3D except 2 factors with pointwise output
            if N==2 && ~(iz==1 || iz==3 && strcmp(kind{1},'pw')) || N==3 && (numel(Z)>2 || strcmp(kind{1},'kron')),     continue,   end
            szlist = {sz1};     if iz==1 || iz==3,      szlist = {sz1,[2;1;2]};     end     % vector-valued states, cost
            for szc = szlist
                sz = szc{1};
                if N>1 && sz(1)>1,      continue,   end
                seed = 10000*N + 1000*numel(kind{1}) + 100*iz + 10*sz(1);
                rng(seed);
                label = sprintf('N=%d,%s,Z=[%s],sz=[%s]',N,kind{1},num2str(Z),num2str(sz'));
                f = make_f(N,Z,sz,kind{1});
                x = struct('a',rand_pol(svs,sz(1)),'b',rand_pol(svs,sz(2)),'c',rand_pol(svs,sz(3)));
                names = allnames(unique(Z));
                % the first state, and all of them, each by a constant of a different type
                nsub = numel(names);
                masks = unique([1,2^nsub-1]);       if N==3,    masks = 1;      end     % cost
                for mask = masks
                    sub = names(logical(bitget(mask,1:numel(names))));
                    lbl = sprintf('%s,sub={%s}',label,strjoin(sub,','));
                    is_all = numel(sub)==numel(names);
                    forms = {'double','polynomial','dpvar'};
                    for form = forms(1+mod(mask+iz,3))
                        cvals = cellfun(@(nm) rand_const(svs,sz(strcmp(allnames,nm)),form{1}),sub,'UniformOutput',false);
                        if is_all && strcmp(kind{1},'kron')
                            expect_error(@() subs(f,sub,vertcat(cvals{:})),'tensor product');
                            ncases = ncases+1;
                        else
                            ncases = ncases+check_const(f,sub,cvals,x,sz,tol,[lbl,',',form{1}]);
                        end
                    end
                    if N<3,     ncases = ncases+check_poly(f,sub,x,sz,N,1+mod(mask,2),tol,lbl);     end
                end
                % each factor by position, by a function of a new state, and by a constant
                poslist = unique([1,numel(Z)]);     if N>1,     poslist = 1;    end         % cost
                if N==3,    poslist = [];   end
                for pos = poslist
                    ncases = ncases+check_pos(f,Z,pos,x,sz,N,tol,[label,sprintf(',pos=%d',pos)]);
                end
            end
        end
    end
end

%% (5b) several monomials and several rows
for N = 1:2
    [svars,dom] = spatial_space(N);     svs = svar_all(1:N);
    rng(20000+N);
    sz = sz1;
    x = struct('a',rand_pol(svs,1),'b',rand_pol(svs,1),'c',rand_pol(svs,1));
    fs = {  make_f(N,[1,2],sz,'pw') + make_f(N,1,sz,'pw') + make_f(N,[2,3],sz,'pw') + make_f(N,[1,2,3],sz,'pw'),   'sum of monomials'
            [make_f(N,[1,2],sz,'pw'); make_f(N,[1,3],sz,'pw')],                                                  'rows'
            [make_f(N,[1,2],sz,'pw') + make_f(N,3,sz,'pw'); make_f(N,[1,1,3],sz,'pw')],                         'rows and monomials' };
    for k=1:size(fs,1)
        f = fs{k,1};
        for sub = {{'a'},{'b'},{'c'},{'a','c'},{'a','b','c'}}
            lbl = sprintf('(5b) N=%d,%s,sub={%s}',N,fs{k,2},strjoin(sub{1},','));
            for form = {'double','polynomial'}
                cvals = cellfun(@(nm) rand_const(svs,1,form{1}),sub{1},'UniformOutput',false);
                ncases = ncases+check_const(f,sub{1},cvals,x,sz,tol,[lbl,',',form{1}]);
            end
            ncases = ncases+check_poly(f,sub{1},x,sz,N,2,tol,lbl);
        end
    end
end

%% (6) functionals
for N = 1:2
    [svars,dom] = spatial_space(N);     svs = svar_all(1:N);
    Zlist_V = {1,[1;2],[1,2],[1,2;3,0],[1,2,0;3,1,2]};
    for iz = 1:numel(Zlist_V)
        Z = Zlist_V{iz};
        if size(Z,2)==3 && N==2,    continue,   end
        rng(30000+10*N+iz);
        V = make_V(svars,dom,Z,0);
        lbl = sprintf('(6) N=%d,Z=%d',N,iz);
        % all states by functions of y
        ncases = ncases+check_poly_V(V,N,tol,lbl);
        % all states by constants: only for monomials of degree one
        cvals = {rand_const(svs,1,'polynomial');rand_const(svs,1,'polynomial');rand_const(svs,1,'polynomial')};
        if size(Z,2)==1
            fnew = subs(V,allnames,vertcat(cvals{:}));
            Yref = eval_polyopvar(V,cvals);
            assert_close(eval_new(fnew,struct(),sz1,allnames),Yref,tol,[lbl,' (all states by constants)']);
            ncases = ncases+1;
        else
            expect_error(@() subs(V,allnames,vertcat(cvals{:})),'degree one');
            ncases = ncases+1;
        end
    end
end

fprintf('Test_polyopvar_subs_class passed (%d checks).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [svars,dom] = spatial_space(N)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
end

function p = rand_pol(svs,m)
% Random m x 1 polynomial vector in the variables svs (cell of pvar)
p = polynomial(randn(m,1));
for i=1:numel(svs)
    p = p + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
end
end

function v = rand_const(svs,m,form)
% Random m x 1 constant function: 'double', 'polynomial' in svs, or 'dpvar' (no decision variables)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = rand_pol(svs,m);
    if strcmp(form,'dpvar'),    v = dpvar(v);   end
end
end

function f = make_f(N,Zrow,sz,kind,raw)
% Function f = Kop*(x_{Z1} o ... o x_{Zd}) in the states a,b,c of sizes sz, all depending on the
% variables s1,...,sN, with Kop random, with pointwise ('pw', output size 2) or Kronecker ('kron') output
[svars,dom] = spatial_space(N);
vars = struct('in',{svars},'out',{svars});
d = numel(Zrow);
din = reshape(sz(Zrow),[],1);
if strcmp(kind,'pw')
    dims = [2*ones(d,1),din];       ot = [1,2];     mout = 2;
else
    dims = [(1+(1:d)'),din];        ot = [2,2];     mout = prod(dims(:,1));
end
dp = struct('in',true(d,N),'out',true(d,N));
Kop = rand_tensopvar(dims,vars,dom,dp,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = Zrow;                 f.varname = {'a';'b';'c'};
f.dims.in = sz(:);          f.dims.out = mout;
f.vars.in = svars;          f.dom.in = dom.in;          f.depmat.in = true(3,N);
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
if nargin<5 || ~raw,      f = initialize(f);      end
end

function V = make_V(svars,dom,Z,q)
% Functional V(a,b,c) = sum_k Aop_k*Z_k, Aop_k random intop, for the monomials in the rows of Z (padded with zeros)
N = numel(svars);
C = [];
for k=1:size(Z,1)
    idcs = Z(k,Z(k,:)~=0);
    Aop = rand_intop([1,1],svars,dom.in,numel(idcs)*ones(1,N),2,q);
    Ck = tensopmat(Aop);
    if k==1,    C = Ck;     else,   C = [C,Ck];     end %#ok<AGROW>
end
V = polyopvar();
V.Z = Z;                    V.C = C;                    V.varname = {'a';'b';'c'};
V.dims.in = ones(3,1);      V.dims.out = 1;
V.vars.in = svars;          V.vars.out = cell(1,0);
V.dom.in = dom.in;          V.dom.out = zeros(0,2);
V.depmat.in = true(3,N);    V.depmat.out = zeros(1,0);
V = initialize(V);
end

function val = make_val(svars,dom,yname,sza,szy,dval)
% val(y)(s): vector of size sza, dval degree-one monomials in the new state y of size szy, random operators
N = numel(svars);
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([sza,szy],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok<AGROW>
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {yname};
val.dims.in = szy;          val.dims.out = sza;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = combine_terms(initialize(val));       % one monomial: vertcat (common_basis) drops duplicate monomials
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

function n = check_const(f,sub,cvals,x,sz,tol,label)
% Substitute the states sub of f by the constants cvals, by name, and by position if f has a single monomial
% in which each is a single factor, and compare with f evaluated at the constants
allnames = {'a','b','c'};
p = randperm(numel(sub));       sub = sub(p);       cvals = cvals(p);       % the order of var and val need not be that of the states
cmat = vertcat(cvals{:});
xs = x;
for k=1:numel(sub)
    xs.(sub{k}) = cvals{k};
    if isa(xs.(sub{k}),'dpvar'),    xs.(sub{k}) = dpvar2poly(xs.(sub{k}));      end
end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
fnew = subs(f,sub,cmat);
assert_close(eval_new(fnew,rmfield(x,sub),sz,allnames),Yref,tol,label+" (by name)");
n = 1;
Zr = full(f.Z);
idx = cellfun(@(s) find(strcmp(allnames,s)),sub);
if size(Zr,1)==1 && all(arrayfun(@(i) sum(Zr==i)==1,idx))
    pos = arrayfun(@(i) find(Zr==i),idx);
    fnew = subs(f,pos(:)',cmat);
    assert_close(eval_new(fnew,rmfield(x,sub),sz,allnames),Yref,tol,label+" (by position)");
    n = 2;
end
end

function svs = pv_cell(N)
% Cell with the pvar objects s1,...,sN
pvar s1 s2 s3
svs = {s1,s2,s3};
svs = svs(1:N);
end

function n = check_poly(f,sub,x,sz,N,dval,tol,label)
% Substitute the states sub of f by functions val_k(y) of a new state y, by name, and compare with f evaluated
% at the values val_k(y)
allnames = {'a','b','c'};
sub = sub(randperm(numel(sub)));        % the order of var and val need not be that of the states
[svars,dom] = spatial_space(N);
svs = pv_cell(N);
szy = 1+mod(N,2);
y = rand_pol(svs,szy);
val = [];       xs = x;
for k=1:numel(sub)
    sza = sz(strcmp(allnames,sub{k}));
    vk = make_val(svars,dom,'b2',sza,szy,dval);
    xs.(sub{k}) = eval_polyopvar(vk,{y});
    if k==1,    val = vk;   else,   val = [val;vk];     end %#ok<AGROW>
end
Yref = eval_polyopvar(f,{xs.a;xs.b;xs.c});
fnew = subs(f,sub,val);
xr = rmfield(x,sub);        xr.b2 = y;
assert_close(eval_new(fnew,xr,[sz;szy],[allnames,{'b2'}]),Yref,tol,label+" (val(y) by name)");
n = 1;
end

function n = check_pos(f,Z,pos,x,sz,N,tol,label)
% Substitute the factor pos of the monomial Z of f (single monomial) by a function of a new state y, and by a
% constant, and compare with Kop applied to the factors in which factor pos is replaced
allnames = {'a','b','c'};
[svars,dom] = spatial_space(N);
svs = pv_cell(N);
szy = 1+mod(N,2);
y = rand_pol(svs,szy);
nm = allnames{Z(pos)};      sza = sz(Z(pos));
fac = arrayfun(@(i) x.(allnames{i}),Z,'UniformOutput',false);
% function of y
val = make_val(svars,dom,'b2',sza,szy,1+mod(pos,2));
fac_y = fac;        fac_y{pos} = eval_polyopvar(val,{y});
Yref = apply_op(f.C,{fac_y});
fnew = subs(f,pos,val);
xr = x;     xr.b2 = y;
assert_close(eval_new(fnew,xr,[sz;szy],[allnames,{'b2'}]),Yref,tol,label+" (val(y))");
% constant
cval = rand_const(svs,sza,'polynomial');
fac_c = fac;        fac_c{pos} = cval;
Yref = apply_op(f.C,{fac_c});
fnew = subs(f,pos,cval);
assert_close(eval_new(fnew,x,sz,allnames),Yref,tol,label+" (constant)");
n = 2;
end

function n = check_poly_V(V,N,tol,label)
% Functional V: all states (by name) by functions of y, evaluated through the kernels (eval_polyopvar)
allnames = {'a','b','c'};
[svars,dom] = spatial_space(N);
svs = pv_cell(N);
y = rand_pol(svs,1);
val = [];       xs = struct();
for k=1:3
    vk = make_val(svars,dom,'b2',1,1,1+mod(k,2));
    xs.(allnames{k}) = eval_polyopvar(vk,{y});
    if k==1,    val = vk;   else,   val = [val;vk];     end %#ok<AGROW>
end
Yref = eval_functional(V,xs);
fnew = subs(V,allnames,val);
Ynew = eval_functional(fnew,struct('b2',y));
assert_close(Ynew,Yref,tol,label+" (all states by val(y))");
n = 1;
end

function Y = eval_functional(V,vals)
% Evaluate the functional V (polyopvar of intop coefficients) at the states vals.(name), polynomials in the spatial
% variables: the factors of each monomial are multiplied (kron), with the variable s_k of the factors that depend
% on it renamed s_k_j, j numbering these factors, as apply_op of 'intop' requires
vnames = V.vars.in;
dep_state = double(V.depmat.in)>0;
Zf = full(V.Z);
Y = 0;
for k=1:size(Zf,1)
    op = V.C.ops{1,k};
    if isempty(op),     continue,   end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    slot = zeros(1,numel(vnames));
    v = 1;
    for j=1:d
        xj = vals.(V.varname{idcs(j)});
        for mm=1:numel(vnames)
            if dep_state(idcs(j),mm)
                slot(mm) = slot(mm)+1;
                if d>1
                    xj.varname(strcmp(xj.varname,vnames{mm})) = {[vnames{mm},'_',num2str(slot(mm))]};
                end
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
end

function assert_close(Ynew,Yref,tol,label)
Yref = polynomial(Yref);        Ynew = polynomial(Ynew);
Yc = Yref.coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_subs_class (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol*max(1,max(abs(full(Yc(:)))))
    error("Test_polyopvar_subs_class failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if isempty(str) || ~isempty(regexp(ME.message,str,'once'))
        return
    end
    error("Test_polyopvar_subs_class failed: unexpected error message '%s', expected '%s'.",ME.message,str)
end
error("Test_polyopvar_subs_class failed: expected an error (%s), but none was thrown.",str)
end
