% Test_polyopvar_subs_multirow.m -- substitution (all routes of
% @polyopvar/subs.m) in polynomials with multiple rows of coefficient
% operators, built by vertcat of two polynomials (f = [f1;f2]), possibly
% with vector-valued states and rows depending on different output
% variables:
%   (1) subs_var by name, by a polynomial function val(y);
%   (2) subs_factor by position, same;
%   (3) subs_const by name and by position, 'double'/'polynomial' vectors.
%
% Verification is independent of subs: the substituted polynomial is
% evaluated at the remaining states and compared with the original f
% evaluated (eval_polyopvar) at a = val(y) or a = constant.
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
sizes = [1,1,1; 2,3,2];             % [ma,mb,my]
seed = 0;

for Ntot = 1:2
    svs = svar_all(1:Ntot);
    % Row dependencies on the output variables: all, or (for Ntot=2) row 1
    % on s1 only and row 2 on s2 only
    outsels = {true(2,Ntot)};
    if Ntot==2,     outsels{2} = [true,false; false,true];      end
    for iout = 1:numel(outsels)
        outsel = outsels{iout};
        for isz = 1:size(sizes,1)
            ma = sizes(isz,1);      mb = sizes(isz,2);      my = sizes(isz,3);
            sz = [ma,mb,my];
            for kind = {'pw','kron','functional'}
                if strcmp(kind{1},'functional') && iout>1,  continue,   end
                for dval = 1:2
                    seed = seed+1;      rng(11000+seed);
                    lbl = sprintf('Ntot=%d,out=%d,[ma,mb,my]=[%d,%d,%d],%s,dval=%d',Ntot,iout,ma,mb,my,kind{1},dval);
                    if strcmp(kind{1},'functional'),    mono = 'a';     else,   mono = 'ab';    end
                    f = vertcat(build_f(Ntot,[ma,mb],kind{1},mono,outsel(1,:)),...
                                build_f(Ntot,[ma,mb],kind{1},mono,outsel(2,:)));
                    if size(f.C.ops,1)~=2
                        error("Test_polyopvar_subs_multirow: test setup failed (%s): f should have 2 rows.",lbl)
                    end
                    val = build_val(Ntot,ma,my,dval);
                    y = rand_pol(svs,my);       b = rand_pol(svs,mb);
                    vy = eval_polyopvar(val,{y});
                    Yref = eval_polyopvar(f,{vy;b});
                    % (1) by name
                    fnew = subs(f,{'a'},val);
                    check(fnew,Yref,{'b',b,'y',y},sz,tol,lbl+" (1)");
                    % (2) by position (first factor)
                    fnew = subs(f,1,val);
                    check(fnew,Yref,{'b',b,'y',y},sz,tol,lbl+" (2)");
                    ncases = ncases+2;
                end
            end

            %% (3) constants
            for cform = {'double','polynomial'}
                rng(12000+100*Ntot+10*iout+isz+numel(cform{1}));
                ca = rand_val(svs,ma,cform{1});     cb = rand_val(svs,mb,cform{1});
                b = rand_pol(svs,mb);
                lbl = sprintf('Ntot=%d,out=%d,[ma,mb]=[%d,%d],const %s',Ntot,iout,ma,mb,cform{1});
                % two-factor pointwise function, all factors replaced
                f = vertcat(build_f(Ntot,[ma,mb],'pw','ab',outsel(1,:)),build_f(Ntot,[ma,mb],'pw','ab',outsel(2,:)));
                Yref = eval_polyopvar(f,{ca;cb});
                fnew = subs(f,{'a';'b'},[ca;cb]);
                check(fnew,Yref,{},sz,tol,lbl+" (3a) by name");
                fnew = subs(f,[1,2],[ca;cb]);
                check(fnew,Yref,{},sz,tol,lbl+" (3a) by position");
                % single-factor function (pointwise, Kronecker) and functional
                for kind = {'pw','kron','functional'}
                    if strcmp(kind{1},'functional') && iout>1,  continue,   end
                    f = vertcat(build_f(Ntot,[ma,mb],kind{1},'a',outsel(1,:)),build_f(Ntot,[ma,mb],kind{1},'a',outsel(2,:)));
                    Yref = eval_polyopvar(f,{ca;b});
                    if strcmp(kind{1},'kron')
                        Yref.varname = regexprep(Yref.varname,"_1$","");
                    end
                    fnew = subs(f,{'a'},ca);
                    check(fnew,Yref,{},sz,tol,lbl+" (3b) "+kind{1});
                    ncases = ncases+1;
                end
                % several monomials, of which only a is substituted: f = [a-part; b-part] terms
                fa = vertcat(build_f(Ntot,[ma,mb],'pw','a',outsel(1,:)),build_f(Ntot,[ma,mb],'pw','a',outsel(2,:)));
                fb = vertcat(build_f(Ntot,[ma,mb],'pw','b',outsel(1,:)),build_f(Ntot,[ma,mb],'pw','b',outsel(2,:)));
                f = fa + fb;
                Yref = eval_polyopvar(f,{ca;b});
                fnew = subs(f,{'a'},ca);
                check(fnew,Yref,{'b',b},sz,tol,lbl+" (3c) kept monomials");
                ncases = ncases+3;
            end
        end
    end
end

fprintf('Test_polyopvar_subs_multirow passed (%d substitution cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function f = build_f(N,sz,kind,mono,outsel)
% Single-row polynomial in vector states a (size sz(1)) and b (size sz(2)),
% both depending on all N spatial variables, with output depending on the
% variables selected by the 1 x N logical OUTSEL:
%   mono 'ab': f = Kop*(a o b);  'a': f = Kop*a;  'b': f = Kop*b
% kind 'pw'/'kron': function with pointwise/Kronecker output;
% kind 'functional' (only 'a'): f = Aop*a, scalar.
ma = sz(1);     mb = sz(2);
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
switch mono
    case 'ab',  Z = [1,2];      din = [ma,mb];
    case 'a',   Z = 1;          din = ma;
    case 'b',   Z = 2;          din = mb;
end
d = numel(Z);
f = polyopvar();
f.Z = Z;                    f.varname = {'a';'b'};
f.dims.in = [ma;mb];
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = true(2,N);
if strcmp(kind,'functional')
    Aop = rand_intop([1,ma],svars,dom.in,ones(1,N),3);
    f.C = tensopmat(Aop);
    f.dims.out = 1;
    f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
else
    outsel = logical(outsel(:)');
    vars = struct('in',{svars},'out',{svars(outsel)});
    domk = struct('in',dom.in,'out',dom.out(outsel,:));
    dep = struct('in',true(d,N),'out',true(d,sum(outsel)));
    if strcmp(kind,'pw')
        dims = [2*ones(d,1),din(:)];        ot = [1,2];     mout = 2;
    else
        dims = [(1+(1:d)'),din(:)];         ot = [2,2];     mout = prod(dims(:,1));
    end
    Kop = rand_tensopvar(dims,vars,domk,dep,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    f.C = tensopmat(Kop);
    f.dims.out = mout;
    f.vars.out = svars(outsel);     f.dom.out = dom.out(outsel,:);      f.depmat.out = Kop.deps.out;
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
    error("Test_polyopvar_subs_multirow (%s): reference value is zero; test would be vacuous.",label)
end
if ~isequal(size(Ynew),size(Yref))
    error("Test_polyopvar_subs_multirow failed (%s): output size [%s] differs from reference [%s].",label,num2str(size(Ynew)),num2str(size(Yref)))
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_subs_multirow failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end
