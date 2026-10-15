% Test_polyopvar_subs_vector_states.m -- substitution (all routes of
% @polyopvar/subs.m) in polynomials of vector-valued states:
%   (1) subs_var by name, by a polynomial function val(y) of a vector
%       state y;
%   (2) subs_factor by position, same;
%   (3) subs_const by name and by position, with 'double'/'polynomial'
%       vectors.
% States a,b,y have sizes ma,mb,my; f is a function (pointwise or
% Kronecker output type) or a functional (degree one), in 1 or 2 spatial
% variables.
%
% Verification is independent of subs: the substituted polynomial is
% evaluated at the remaining states and compared with the original f
% evaluated (eval_polyopvar, apply_op) at a = val(y) or a = constant.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.

pvar s1 s2
svar_all = {s1,s2};
tol = 1e-8;
ncases = 0;
sizes = [2,1,1; 1,3,2; 2,3,2];          % [ma,mb,my]
seed = 0;

for Ntot = 1:2
    svs = svar_all(1:Ntot);
    for isz = 1:size(sizes,1)
        ma = sizes(isz,1);      mb = sizes(isz,2);      my = sizes(isz,3);
        for kind = {'pw','kron','functional'}
            for dval = 1:2
                seed = seed+1;      rng(7000+seed);
                lbl = sprintf('Ntot=%d,[ma,mb,my]=[%d,%d,%d],%s,dval=%d',Ntot,ma,mb,my,kind{1},dval);
                if strcmp(kind{1},'functional')
                    f = build_f(Ntot,[ma,mb],kind{1},'a');
                else
                    f = build_f(Ntot,[ma,mb],kind{1},'ab');
                end
                val = build_val(Ntot,ma,my,dval);
                y = rand_pol(svs,my);       b = rand_pol(svs,mb);
                vy = eval_polyopvar(val,{y});
                if strcmp(kind{1},'functional')
                    % degree-one functional in a: f(a,b) = Aop*a
                    f = build_f(Ntot,[ma,mb],kind{1},'a');
                    Yref = eval_polyopvar(f,{vy;b});
                    % (1) by name
                    fnew = subs(f,{'a'},val);
                    check(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (1)");
                    % (2) by position
                    fnew = subs(f,1,val);
                    check(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (2)");
                    ncases = ncases+2;
                else
                    Yref = eval_polyopvar(f,{vy;b});
                    % (1) by name
                    fnew = subs(f,{'a'},val);
                    check(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (1)");
                    % (2) by position (first factor of a o b)
                    fnew = subs(f,1,val);
                    check(fnew,Yref,{'b',b,'y',y},[ma,mb,my],tol,lbl+" (2)");
                    ncases = ncases+2;
                end
            end
        end

        %% (2b) by position in a o a: replace one factor, keep the other
        for kind = {'pw','kron'}
            rng(8000+10*Ntot+isz);
            lbl = sprintf('Ntot=%d,[ma,mb,my]=[%d,%d,%d],%s,a o a',Ntot,ma,mb,my,kind{1});
            [f,Kop] = build_f(Ntot,[ma,mb],kind{1},'aa');
            val = build_val(Ntot,ma,my,1);
            x = rand_pol(svs,ma);       y = rand_pol(svs,my);
            vy = eval_polyopvar(val,{y});
            for fctr = 1:2
                X = {x,x};      X{fctr} = vy;
                Yref = apply_op(Kop,X);
                fnew = subs(f,fctr,val);
                check(fnew,Yref,{'a',x,'y',y},[ma,mb,my],tol,lbl+sprintf(" fctr=%d",fctr));
                ncases = ncases+1;
            end
        end

        %% (3) constants, by name and by position
        for cform = {'double','polynomial'}
            rng(9000+10*Ntot+isz+numel(cform{1}));
            ca = rand_val(svs,ma,cform{1});     cb = rand_val(svs,mb,cform{1});
            b = rand_pol(svs,mb);
            lbl = sprintf('Ntot=%d,[ma,mb]=[%d,%d],const %s',Ntot,ma,mb,cform{1});
            % two-factor pointwise function, all factors replaced
            f = build_f(Ntot,[ma,mb],'pw','ab');
            Yref = eval_polyopvar(f,{ca;cb});
            fnew = subs(f,{'a';'b'},[ca;cb]);
            check(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3a) by name");
            fnew = subs(f,[1,2],[ca;cb]);
            check(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3a) by position");
            fnew = subs(f,{'b';'a'},[cb;ca]);
            check(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3a) reordered names");
            % single-factor function and functional
            for kind = {'pw','kron','functional'}
                f = build_f(Ntot,[ma,mb],kind{1},'a');
                Yref = eval_polyopvar(f,{ca;b});
                if strcmp(kind{1},'kron')
                    Yref.varname = regexprep(Yref.varname,"_1$","");
                end
                fnew = subs(f,{'a'},ca);
                check(fnew,Yref,{},[ma,mb,my],tol,lbl+" (3b) "+kind{1});
                ncases = ncases+1;
            end
            % partial substitution: the factor replaced remains as a known factor
            f = build_f(Ntot,[ma,mb],'pw','ab');
            Yref = eval_polyopvar(f,{ca;b});
            fnew = subs(f,{'a'},ca);
            check(fnew,Yref,{'b',b},[ma,mb,my],tol,lbl+" (3c) partial, by name");
            fnew = subs(f,1,ca);
            check(fnew,Yref,{'b',b},[ma,mb,my],tol,lbl+" (3c) partial, by position");
            ncases = ncases+2;
            % wrong number of rows (scalar for a vector state) errors
            if ma>1
                expect_error(@() subs(f,{'a';'b'},[ca(1);cb]),'Number of rows');
            end
            ncases = ncases+5;
        end
    end
end

fprintf('Test_polyopvar_subs_vector_states passed (%d substitution cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function [f,Kop] = build_f(N,sz,kind,mono)
% Polynomial in vector states a (size sz(1)) and b (size sz(2)), both
% depending on all N spatial variables:
%   mono 'ab': f = Kop*(a o b);   'aa': f = Kop*(a o a);   'a': f = Kop*a
% kind 'pw'/'kron': function with pointwise/Kronecker output (output size 2
% resp. 2*2 for the pointwise/Kronecker product of two factors; 2 for
% one factor); kind 'functional' (only mono 'a'): f = Aop*a, scalar.
ma = sz(1);     mb = sz(2);
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
switch mono
    case 'ab',  Z = [1,2];      din = [ma,mb];
    case 'aa',  Z = [1,1];      din = [ma,ma];
    case 'a',   Z = 1;          din = ma;
end
d = numel(Z);
f = polyopvar();
f.Z = Z;                    f.varname = {'a';'b'};
f.dims.in = [ma;mb];
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = true(2,N);
Kop = [];
if strcmp(kind,'functional')
    Aop = rand_intop([1,ma],svars,dom.in,ones(1,N),3);
    f.C = tensopmat(Aop);
    f.dims.out = 1;
    f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
else
    dep = struct('in',true(d,N),'out',true(d,N));
    if strcmp(kind,'pw')
        dims = [2*ones(d,1),din(:)];        ot = [1,2];     mout = 2;
    else
        dims = [(1+(1:d)'),din(:)];         ot = [2,2];     mout = prod(dims(:,1));
    end
    Kop = rand_tensopvar(dims,vars,dom,dep,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    f.C = tensopmat(Kop);
    f.dims.out = mout;
    f.vars.out = svars;     f.dom.out = dom.out;    f.depmat.out = Kop.deps.out;
end
f = initialize(f);
end

function val = build_val(N,ma,my,dval)
% val(y) in R^ma, dval degree-one monomials in the vector state y (size my)
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
dep1 = struct('in',true(1,N),'out',true(1,N));
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([ma,my],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = my;           val.dims.out = ma;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol(svs,m)
p = polynomial(randn(m,1));
for i=1:numel(svs)
    p = p + randn(m,1)*svs{i} + randn(m,1)*svs{i}^2;
end
end

function v = rand_val(svs,m,form)
if strcmp(form,'double')
    v = randn(m,1);
else
    v = rand_pol(svs,m);
end
end

function check(fnew,Yref,vals,sz,tol,label)
% Evaluate fnew at states given as {'name',value,...}; unused states of
% fnew (e.g. the substituted-out 'a') get zeros of the right size
sizes = struct('a',sz(1),'b',sz(2),'y',sz(3));
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    xv{k} = polynomial(zeros(sizes.(fnew.varname{k}),1));
    for i=1:2:numel(vals)
        if strcmp(fnew.varname{k},vals{i}),     xv{k} = vals{i+1};     end
    end
end
Ynew = eval_polyopvar(fnew,xv);
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_subs_vector_states (%s): reference value is zero; test would be vacuous.",label)
end
if ~isequal(size(Ynew),size(Yref))
    error("Test_polyopvar_subs_vector_states failed (%s): output size [%s] differs from reference [%s].",label,num2str(size(Ynew)),num2str(size(Yref)))
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_subs_vector_states failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_polyopvar_subs_vector_states failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_polyopvar_subs_vector_states failed: expected an error containing '%s', but none was thrown.",str)
end
