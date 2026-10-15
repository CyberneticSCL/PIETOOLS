% Test_polyopvar_subs_var_function.m -- substitution of a state variable
% by name (subs_var in @polyopvar/subs.m) in
%   (1) polynomial functions f(a,b)(s) = Kop*(a o b), tensopvar
%       coefficient,
%   (2) polynomial functionals f(a,b) = Aop*(a o b), intop coefficient,
% by (i) a polynomial function val(y) in a new state y, (ii) a constant.
% The variable is given as cellstr or as 'polyopvar' (polyopvar('a',...)).
%
% Verification is semantic and independent of subs: fnew(b,y) is compared
% to the original f evaluated (eval_polyopvar) at a=val(y) -- itself
% evaluated by eval_polyopvar -- and the unchanged b.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Evaluation uses eval_polyopvar (test-only, does not use subs), not
%   apply_poly, which evaluates by subs and would make the check circular.
% - Kronecker outputs: the reference evaluation renames output variables
%   (suffix "_1"); stripped before comparing in (ii).

pvar s1 s2
svar_all = {s1,s2};
tol = 1e-8;
ncases = 0;

%% (i) substitution by a polynomial function of a new state y
for Ntot = 1:2
    for kind = {'function_pw','function_kron','functional'}
        for dval = 1:2
            for varform = {'cellstr','polyopvar'}
                rng(1000*Ntot+100*dval+10*numel(kind{1})+numel(varform{1}));
                [f,val,svars,dom] = build_f(Ntot,kind{1},false,dval);
                if strcmp(varform{1},'cellstr')
                    var = {'a'};
                else
                    var = polyopvar('a',svars,dom.in);
                end
                fnew = subs(f,var,val);
                check_subs(f,fnew,val,svar_all(1:Ntot),tol,...
                    sprintf('(i) Ntot=%d,%s,dval=%d,var=%s',Ntot,kind{1},dval,varform{1}));
                ncases = ncases+1;
            end
        end
    end
end

%% (ii) substitution by a constant, for a state a not depending on any spatial variable
% The constant replaces the entire (degree-one) monomial, giving a degree-0 polynomial,
% or some factors of a monomial of higher degree, which remain as factors that are multipliers
% from R (see polyopvar), also for a 'double' constant (partial substitution in subs_const).
% Verified against the original polynomial evaluated at the constant.
for Ntot = 1:2
    for ot = {[1,2],[2,2]}
        for cform = {'polyopvar','double'}
            rng(5000+100*Ntot+10*ot{1}(1)+numel(cform{1}));
            c = randn;
            if strcmp(cform{1},'polyopvar')
                val = constant2polyopvar(c);
            else
                val = c;
            end
            % Degree one in a: f(a)(s) = Kop*a
            f = build_const_f(Ntot,1,ot{1});
            fnew = subs(f,{'a'},val);
            Yref = eval_polyopvar(f,{polynomial(c)});
            if strcmp(cform{1},"double")    % subs_const returns the output in the original variables
                Yref.varname = regexprep(Yref.varname,"_1$","");   % apply_tensopvar renames output variables for Kronecker outputs
            end
            Ynew = eval_polyopvar(fnew,cell(numel(fnew.varname),1));
            assert_close(Ynew,Yref,tol,sprintf('(ii) Ntot=%d,ot=[%d,%d],val=%s',Ntot,ot{1}(1),ot{1}(2),cform{1}));
            ncases = ncases+1;
            % Degree two in a: a constant polyopvar replaces both factors, leaving a degree-0
            % term with two factors (a 'double' constant is handled by subs_const instead,
            % see Test_polyopvar_subs_const)
            if strcmp(cform{1},'polyopvar')
                f2 = build_const_f(Ntot,2,ot{1});
                fnew2 = subs(f2,{'a'},val);
                Yref2 = eval_polyopvar(f2,{polynomial(c)});
                Ynew2 = eval_polyopvar(fnew2,cell(numel(fnew2.varname),1));
                assert_close(Ynew2,Yref2,tol,sprintf('(ii) a o a, Ntot=%d,ot=[%d,%d]',Ntot,ot{1}(1),ot{1}(2)));
                if any(any(fnew2.Z))
                    error("Test_polyopvar_subs_var_function failed: a o a with a constant should be a degree-0 term.")
                end
                ncases = ncases+1;
            end
            % Constant in a monomial also containing b: the constant remains as a factor, Z=[0,b],
            % for a constant polyopvar as well as for a 'double' (partial substitution in subs_const)
            for kind3 = {'function_pw','function_kron'}
                f3 = build_f(Ntot,kind3{1},true,1);
                fnew3 = subs(f3,{'a'},val);
                b3 = rand_pol(svar_all(1:Ntot));
                Yref3 = eval_polyopvar(f3,{polynomial(c);b3});
                xv3 = cell(numel(fnew3.varname),1);
                for k3=1:numel(xv3)
                    if strcmp(fnew3.varname{k3},'b'),   xv3{k3} = b3;
                    else,                               xv3{k3} = polynomial(0);    end
                end
                assert_close(eval_polyopvar(fnew3,xv3),Yref3,tol,sprintf('(ii) a o b, Ntot=%d,%s,val=%s',Ntot,kind3{1},cform{1}));
                if ~isequal(full(fnew3.Z),[0,2]) && ~isequal(full(fnew3.Z),[0,1])   % constant, then b
                    error("Test_polyopvar_subs_var_function failed: a o b with a constant should have Z=[0,b], not [%s].",num2str(full(fnew3.Z)))
                end
                ncases = ncases+1;
            end
        end
    end
end

fprintf('Test_polyopvar_subs_var_function passed (%d substitution cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_subs(f,fnew,val,svs,tol,label)
a = rand_pol(svs);      b = rand_pol(svs);      y = rand_pol(svs);
vy = eval_polyopvar(val,{y});
Yref = eval_polyopvar(f,{vy;b});
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    if strcmp(fnew.varname{k},'y'),     xv{k} = y;
    elseif strcmp(fnew.varname{k},'b'), xv{k} = b;
    elseif strcmp(fnew.varname{k},'a')
        % Substituted-out state may remain declared, but must not appear in Z
        if any(fnew.Z(:)==k)
            error("Test_polyopvar_subs_var_function failed (%s): fnew.Z still references the substituted state 'a'.",label)
        end
        xv{k} = polynomial(0);
    else
        error("Test_polyopvar_subs_var_function failed (%s): unexpected state '%s' in result.",label,fnew.varname{k})
    end
end
Ynew = eval_polyopvar(fnew,xv);
assert_close(Ynew,Yref,tol,label)
end

function assert_close(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_subs_var_function (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_subs_var_function failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function [f,val,svars,dom] = build_f(N,kind,a_const,dval)
% f(a,b) with monomial a o b; if a_const, state a does not depend on any
% spatial variable. val(y): dval degree-one monomials in a new state y.
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
depin = [true(1,N)&~a_const; true(1,N)];
f = polyopvar();
f.Z = [1,2];                f.varname = {'a';'b'};
f.dims.in = [1;1];          f.dims.out = 1;
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = depin;
switch kind
    case {'function_pw','function_kron'}
        if strcmp(kind,'function_pw'),  ot = [1,2];     else,   ot = [2,2];     end
        dep = struct('in',depin,'out',true(2,N));
        Kop = rand_tensopvar(ones(2,2),vars,dom,dep,1,ot,2);
        if ~Kop.chkval,     Kop = initialize(Kop);      end
        f.C = tensopmat(Kop);
        f.vars.out = svars;     f.dom.out = dom.out;    f.depmat.out = Kop.deps.out;
    case 'functional'
        % degree-one functional only: f(a)=Aop*a
        f.Z = 1;
        Aop = rand_intop([1,1],svars,dom.in,double(depin(1,:)),3);
        f.C = tensopmat(Aop);
        f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
end
f = initialize(f);
% val
dep1 = struct('in',true(1,N),'out',depin(1,:));
if a_const
    val = [];   return
end
dep1.out = true(1,N);
C = [];
for k=1:dval
    Bk = tensopmat(rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],2));
    if k==1,    C = Bk;     else,   C = [C,Bk];     end %#ok
end
val = polyopvar();
val.Z = ones(dval,1);       val.C = C;                  val.varname = {'y'};
val.dims.in = 1;            val.dims.out = 1;
val.vars.in = svars;        val.vars.out = svars;
val.dom.in = dom.in;        val.dom.out = dom.out;
val.depmat.in = true(1,N);  val.depmat.out = ones(1,N);
val = initialize(val);
end

function p = rand_pol(svs)
p = polynomial(randn);
for i=1:numel(svs)
    p = p + randn*svs{i} + randn*svs{i}^2;
end
end

function f = build_const_f(N,d,ot)
% f(a)(s) = Kop*(a o ... o a) (d factors), a not depending on spatial variables
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',zeros(0,2),'out',repmat([0,1],N,1));
vars = struct('in',{cell(1,0)},'out',{svars});
dep = struct('in',false(d,0),'out',true(d,N));
Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
if ~Kop.chkval,     Kop = initialize(Kop);      end
f = polyopvar();
f.Z = ones(1,d);            f.varname = {'a'};
f.dims.in = 1;              f.dims.out = 1;
f.vars.in = cell(1,0);      f.dom.in = zeros(0,2);      f.depmat.in = false(1,0);
f.C = tensopmat(Kop);
f.vars.out = svars;         f.dom.out = dom.out;        f.depmat.out = Kop.deps.out;
f = initialize(f);
end
