% Test_polyopvar_subs_name.m -- substitution of states by name, and of
% several factors at once, in polynomial functionals
% f = sum_k Aop_k*Z_k(x) (intop coefficients, as V(x,u) of piesos_LF),
% see subs_var, subs_fctrs, subs_factor in @polyopvar/subs.m:
%   (i)-(iii) the states a and/or b of f(a,b) replaced by a function
%         val(y) of a new state y, by name, for monomials with repeated
%         factors, e.g. a o a, and several states, a o b o a (i), of
%         different degree, with Z padded with zeros (ii), and kernels
%         that are 'dpvar' with a different list of decision variables
%         in each region of integration, as produced by piesos_quadvar
%         (iii);
%   (iv)  several factors in a monomial, by number;
%   (v)   the new value containing the substituted state itself,
%         a --> val(a);
%   (vi)  renaming of states, also to the name of another state;
%   (vii) repeated factors in polynomial functions (tensopvar
%         coefficients), f(x) = Kop*(x o x);
%   errors for repeated or too large factor numbers, and for states on
%   different function spaces.
%
% Verification is semantic and independent of subs: fnew(y,b) is compared
% to the original f evaluated at the states a=val(y), b (evaluating the
% functional by integration, apply_op, of the kernels at the products of
% the polynomial states), for random polynomial states.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - All kernels of an 'intop' share one list of decision variables Zd, so
%   "a different list in each region" is realized by thin_dvars: a random
%   subset of Zd with nonzero coefficients.
% - States are evaluated with eval_polyopvar (Testfolder/polyopvar), not
%   apply_poly, since apply_poly evaluates by subs, the routine under
%   test.
% - eval_functional must name the dummy variables s_k_j exactly as
%   apply_op does (matched by name).
% - assert_close and (vii) reject a near-zero reference value, so no case
%   can pass vacuously.

pvar s1 s2
svar_all = {s1,s2};
tol = 1e-8;
ncases = 0;

%% (i)-(iii) substitution of one or both states of f(a,b), for a function val(y), by name
%   (Z rows padded with zeros mix the degrees, rows with repeated a, b)
Zlist = {   [1,1];
            [1,2;  2,1];
            [1,1,2;  2,0,0;  1,2,0];
            [1,1,2;  2,2,1;  1,1,1]   };
for Ntot = 1:2
    for iz = 1:numel(Zlist)
        Z = Zlist{iz};
        if Ntot==2 && size(Z,2)>2,  continue,   end         % skip monomials of 3 factors in 2D
        for q = [0,7]
            for dval = 1:2
                seed = 1000*Ntot + 100*iz + 10*dval + q;
                rng(seed);
                [svars,dom] = spatial_space(Ntot);
                f = build_f(svars,dom,{'a';'b'},Z,2,q);
                val = [build_val(svars,dom,'y',dval); build_val(svars,dom,'y',dval)];
                svs = svar_all(1:Ntot);
                d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
                label = sprintf('Ntot=%d,Z=%d,q=%d,dval=%d',Ntot,iz,q,dval);
                y = rand_pol(svs);      b = rand_pol(svs);
                va = eval_polyopvar(val(1,:),{y});     vb = eval_polyopvar(val(2,:),{y});

                % (a) substitute a only: val is a single row
                fnew = subs(f,{'a'},val(1,:));
                no_state_left(fnew,'a',label);
                Yref = eval_functional(f,struct('a',va,'b',b),d0map);
                Ynew = eval_functional(fnew,struct('y',y,'b',b),d0map);
                assert_close(Ynew,Yref,tol,['(a) ',label]);
                ncases = ncases+1;

                % (b) substitute both a and b, listed in a different order than the states
                fnew = subs(f,{'b','a'},[val(2,:);val(1,:)]);
                no_state_left(fnew,'a',label);       no_state_left(fnew,'b',label);
                Yref = eval_functional(f,struct('a',va,'b',vb),d0map);
                Ynew = eval_functional(fnew,struct('y',y),d0map);
                assert_close(Ynew,Yref,tol,['(b) ',label]);
                ncases = ncases+1;
            end
        end
    end
end

%% (iv) substitution of several factors in a single monomial by index
for Ntot = 1:2
    for Zc = {[1,2,0], [1,1], [2,1,2]}
        Z = Zc{1};
        if Ntot==2 && size(Z,2)>2,  continue,   end
        for q = [0,5]
            rng(2000+10*Ntot+q+size(Z,2));
            [svars,dom] = spatial_space(Ntot);
            svs = svar_all(1:Ntot);
            f = build_f(svars,dom,{'a';'b'},Z,2,q);
            d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
            y = rand_pol(svs);      b = rand_pol(svs);
            nfac = nnz(Z);
            val = build_val(svars,dom,'y',1);
            for k=2:nfac,   val = [val; build_val(svars,dom,'y',1)];     end %#ok<AGROW>
            vy = cell(nfac,1);
            for k=1:nfac,   vy{k} = eval_polyopvar(val(k,:),{y});    end
            label = sprintf('Ntot=%d,Z=[%s],q=%d',Ntot,num2str(Z),q);
            % all factors
            fnew = subs(f,1:nfac,val);
            Yref = eval_functional_fctrs(f,Z,vy,1:nfac,struct('a',b,'b',b),d0map);
            Ynew = eval_functional(fnew,struct('y',y),d0map);
            assert_close(Ynew,Yref,tol,['(iv) all factors, ',label]);
            ncases = ncases+1;
            % the last factor only, the others remain states a, b (both given the value b)
            fnew = subs(f,nfac,val(nfac,:));
            Yref = eval_functional_fctrs(f,Z,vy(nfac),nfac,struct('a',b,'b',b),d0map);
            Ynew = eval_functional(fnew,struct('y',y,'a',b,'b',b),d0map);
            assert_close(Ynew,Yref,tol,['(iv) last factor, ',label]);
            ncases = ncases+1;
        end
    end
end

%% (v) the new value contains the substituted state itself, a --> val(a), so that all of the
%   factors must be substituted at once (the substituted factors may not be substituted again)
for q = [0,5]
    for Zc = {[1,1], [1,2,1]}
        Z = Zc{1};
        rng(2500+q+size(Z,2));
        [svars,dom] = spatial_space(1);
        svs = svar_all(1);
        f = build_f(svars,dom,{'a';'b'},Z,2,q);
        d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
        a = rand_pol(svs);      b = rand_pol(svs);
        val = build_val(svars,dom,'a',2);
        va = eval_polyopvar(val,{a});
        label = sprintf('Z=[%s],q=%d',num2str(Z),q);
        fnew = subs(f,{'a'},val);
        Yref = eval_functional(f,struct('a',va,'b',b),d0map);
        Ynew = eval_functional(fnew,struct('a',a,'b',b),d0map);
        assert_close(Ynew,Yref,tol,['(v) a --> val(a), ',label]);
        ncases = ncases+1;
    end
end

%% (vi) renaming of states. For Ntot=2, the state a depends on s1 and s2, but b only on s1,
%   so that the function spaces must move with the names
for Ntot = 1:2
    for q = [0,5]
        rng(3000+10*Ntot+q);
        [svars,dom] = spatial_space(Ntot);
        svs = svar_all(1:Ntot);
        depmat = true(2,Ntot);
        if Ntot==2,     depmat(2,2) = false;    end
        Zren = [1,2,0;  1,1,0;  1,1,2];         % each monomial contains a, so depends on all of the variables
        f = build_f(svars,dom,{'a';'b'},Zren(1:(2+(Ntot==1)),:),2,q,depmat);
        d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:max(q,1),'UniformOutput',false),num2cell(randn(max(q,1),1)));
        a = rand_pol(svs);      b = rand_pol(svs(depmat(2,:)));
        label = sprintf('Ntot=%d,q=%d',Ntot,q);
        % a --> z moves the state after b
        fnew = subs(f,{'a'},{'z'});
        assert(isequal(fnew.varname,{'b';'z'}),"Test_polyopvar_subs_name (%s): states renamed or not sorted.",label)
        Yref = eval_functional(f,struct('a',a,'b',b),d0map);
        Ynew = eval_functional(fnew,struct('z',a,'b',b),d0map);
        assert_close(Ynew,Yref,tol,['(vi) a --> z, ',label]);
        ncases = ncases+1;
        % a,b --> b,a swaps the states
        fnew = subs(f,{'a','b'},{'b','a'});
        Yref = eval_functional(f,struct('a',a,'b',b),d0map);
        Ynew = eval_functional(fnew,struct('b',a,'a',b),d0map);
        assert_close(Ynew,Yref,tol,['(vi) a,b --> b,a, ',label]);
        ncases = ncases+1;
        if Ntot==1
            % a --> b merges the states
            fnew = subs(f,{'a'},{'b'});
            assert(isequal(fnew.varname,{'b'}),"Test_polyopvar_subs_name (%s): states not merged.",label)
            Yref = eval_functional(f,struct('a',b,'b',b),d0map);
            Ynew = eval_functional(fnew,struct('b',b),d0map);
            assert_close(Ynew,Yref,tol,['(vi) a --> b, ',label]);
            ncases = ncases+1;
        else
            % states on different function spaces cannot be merged
            expect_error(@() subs(f,{'a'},{'b'}),'same function space');
        end
    end
end

%% (vii) polynomial FUNCTIONS (tensopvar coefficients) with repeated factors, f(x) = Kop*(x o ... o x),
%   x --> val(y) by name, for pointwise and Kronecker outputs, compared to Kop applied to val(y)
for Ntot = 1:2
    for d = 2:(4-Ntot)          % evaluating the products of 3 factors in 2D takes long
        for otc = {[1,2],[2,2]}
            ot = otc{1};
            for dval = 1:2
                rng(100*Ntot+10*d+ot(1)+dval);
                [f,val,Kop] = build_function_case(Ntot,d,ot,dval);
                x = rand_pol(svar_all(1:Ntot));         y = rand_pol(svar_all(1:Ntot));
                vy = eval_polyopvar(val,{y});
                Yref = apply_op(Kop,repmat({vy},1,d));
                fnew = subs(f,{'x'},val);
                xv = cell(numel(fnew.varname),1);
                for k=1:numel(xv)
                    if strcmp(fnew.varname{k},'y'),     xv{k} = y;
                    else,                               xv{k} = x;      end
                end
                Ynew = eval_polyopvar(fnew,xv);
                label = sprintf('(vii) function, Ntot=%d,d=%d,otype=[%d,%d],dval=%d',Ntot,d,ot(1),ot(2),dval);
                Yc = polynomial(Yref).coeff;
                if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
                    error("Test_polyopvar_subs_name (%s): reference value is zero; test would be vacuous.",label)
                end
                err = polynomial(Ynew-Yref);
                if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
                    error("Test_polyopvar_subs_name failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
                end
                ncases = ncases+1;
            end
        end
    end
end

%% Errors
rng(4000);
[svars,dom] = spatial_space(1);
f = build_f(svars,dom,{'a';'b'},[1,2],2,0);
val = build_val(svars,dom,'y',1);
expect_error(@() subs(f,[1,1],[val;val]),'unique');
expect_error(@() subs(f,[1,3],[val;val]),'degree');
% states on different function spaces cannot be merged; depmat declared consistently with the
% operators (assigning depmat.in alone sets chkval=false and initialize rejects the object)
fc = build_f(svars,dom,{'a';'b'},[1,2],2,0,[true;false]);
expect_error(@() subs(fc,{'a'},{'b'}),'same function space');

fprintf('Test_polyopvar_subs_name passed (%d substitution cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [svars,dom] = spatial_space(N)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
end

function f = build_f(svars,dom,names,Z,pdegs,q,depmat)
% Functional f(a,b) = sum_k Aop_k*Z_k, for monomials in the rows of Z (padded with zeros),
% states depending on the spatial variables indicated by depmat (default: all). If q>0, the
% kernels are 'dpvar' in q decision variables, with a different, random subset of them in
% each region of integration
N = numel(svars);
nx = numel(names);
if nargin<7
    depmat = true(nx,N);
end
C = [];
for k=1:size(Z,1)
    idcs = Z(k,Z(k,:)~=0);
    Aop = rand_intop([1,1],svars,dom.in,sum(depmat(idcs,:),1),pdegs,q);
    if q>0
        Aop = thin_dvars(Aop);
    end
    Ck = tensopmat(Aop);
    if k==1,    C = Ck;     else,   C = [C,Ck];     end %#ok<AGROW>
end
f = polyopvar();
f.Z = Z;                    f.C = C;                    f.varname = names;
f.dims.in = ones(nx,1);     f.dims.out = 1;
f.vars.in = svars;          f.vars.out = cell(1,0);
f.dom.in = dom.in;          f.dom.out = zeros(0,2);
f.depmat.in = depmat;        f.depmat.out = zeros(1,0);
f = initialize(f);
end

function Aop = thin_dvars(Aop)
% Retain a different random subset of the decision variables in each parameter of Aop
% (the rows of the others, in the blocks of rows of each row of the matrix, are set to zero)
q = numel(Aop.Zd);
m = size(Aop,1);
for k=1:numel(Aop.params)
    P = Aop.params{k};
    if q==0,        continue,       end
    J = sort(randperm(q,ceil(q/2)));
    rows = reshape(((0:m-1)'*(q+1) + [1,1+J]).',[],1);
    P(setdiff(1:size(P,1),rows),:) = 0;
    Aop.params{k} = P;
end
end

function val = build_val(svars,dom,yname,dval)
% val(y)(s): dval degree-one monomials in the new state y, random operators
N = numel(svars);
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok<AGROW>
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {yname};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function [f,val,Kop] = build_function_case(N,d,ot,dval)
% Function f(x)(s) = Kop*(x o ... o x), d factors, and val(y) with dval degree-one monomials
% in the new state y (as Test_polyopvar_subs_fctr_function.m)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep = struct('in',true(d,N),'out',true(d,N));
Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = ones(1,d);            f.C = tensopmat(Kop);       f.varname = {'x'};
f.dims.in = 1;              f.dims.out = 1;
f.vars.in = svars;          f.vars.out = svars;
f.dom.in = dom.in;          f.dom.out = dom.out;
f.depmat.in = true(1,N);    f.depmat.out = Kop.deps.out;
f = initialize(f);
val = build_val(svars,dom,'y',dval);
end

function p = rand_pol(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function Y = eval_functional(V,vals,d0map)
% Evaluate the functional V (polyopvar of intop kernels) at the states vals.(name),
% polynomials in the spatial variables. The dummy variables of s_k are numbered over the
% factors of the monomial that depend on s_k, and named s_k_j (as in apply_op)
Zf = full(V.Z);
fac = cell(size(Zf,1),1);
for k=1:size(Zf,1)
    idcs = Zf(k,Zf(k,:)~=0);
    fac{k} = cellfun(@(nm) vals.(nm),reshape(V.varname(idcs),1,[]),'UniformOutput',false);
end
Y = eval_terms(V,fac,d0map);
end

function Y = eval_functional_fctrs(V,Z,vals,fctrs,rest,d0map)
% Evaluate the functional V, defined by the single monomial Z, with the factors fctrs
% replaced by the polynomials vals{:}, and the others by the states in the struct rest
idcs = Z(Z~=0);
fac = cell(1,numel(idcs));
for j=1:numel(idcs)
    fac{j} = rest.(V.varname{idcs(j)});
end
for k=1:numel(fctrs)
    fac{fctrs(k)} = vals{k};
end
Y = eval_terms(V,{fac},d0map);
end

function Y = eval_terms(V,fac,d0map)
% Sum over the monomials k of Aop_k applied to the product of the factors fac{k} (cell of
% polynomials) of the monomial
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
        xj = fac{k}{j};
        for mm=1:numel(vnames)
            if dep_state(idcs(j),mm)
                slot(mm) = slot(mm)+1;
                if d>1
                    isv = strcmp(xj.varname,vnames{mm});
                    xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                end
            end
        end
        v = kron(v,xj);
    end
    Y = Y + apply_op(op,v);
end
if isa(Y,'dpvar')
    Y = eval_dpvar_at(Y,d0map);
end
Y = double(polynomial(Y));
end

function no_state_left(f,name,label)
% The substituted state may remain declared, but must not appear in Z
[~,k] = ismember(name,f.varname);
if k>0 && any(f.Z(:)==k)
    error("Test_polyopvar_subs_name failed (%s): Z still references the substituted state '%s'.",label,name)
end
end

function assert_close(Ynew,Yref,tol,label)
if abs(Yref)<1e-6
    error("Test_polyopvar_subs_name (%s): reference value %g is zero; test would be vacuous.",label,Yref)
end
if abs(Ynew-Yref)>tol*max(1,abs(Yref))
    error("Test_polyopvar_subs_name failed (%s): value %.10g, expected %.10g.",label,Ynew,Yref)
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_polyopvar_subs_name failed: unexpected error message '%s'.",ME.message)
    end
    return
end
error("Test_polyopvar_subs_name failed: expected an error, but none was thrown.")
end
