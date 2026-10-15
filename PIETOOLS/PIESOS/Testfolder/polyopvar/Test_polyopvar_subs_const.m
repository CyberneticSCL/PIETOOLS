% Test_polyopvar_subs_const.m -- substitution by a constant (subs_const
% in @polyopvar/subs.m), i.e. 'double' or 'polynomial' values for
% states/factors, by name (cellstr) and by factor index. If only some of
% the factors in a monomial are substituted, the others remain and the
% constant becomes a known factor (part (6)).
%
% Oracle: fnew is evaluated at the remaining states (eval_polyopvar) and
% compared with the original f evaluated at the constants and the same
% remaining states.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Evaluation uses eval_polyopvar (Testfolder/polyopvar), not apply_poly:
%   apply_poly evaluates by subs, which this test checks.

pvar s1 s2
svar_all = {s1,s2};
tol = 1e-8;
ncases = 0;

for Ntot = 1:2
    svs = svar_all(1:Ntot);
    for cform = {'double','polynomial'}
        rng(300+10*Ntot+numel(cform{1}));
        label = string(sprintf('Ntot=%d,%s',Ntot,cform{1}));
        ca = rand_val(svs,cform{1});       cb = rand_val(svs,cform{1});
        b = rand_val(svs,'polynomial');

        %% (1) All factors of a two-factor function monomial, by name
        f = build_f(Ntot,[1,2],'ab');
        fnew = subs(f,{'a';'b'},[ca;cb]);
        assert_deg0(fnew,label+" (1)");
        assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (1)");
        % ... the order of the names must be respected
        fnew = subs(f,{'b';'a'},[cb;ca]);
        assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (1b)");
        ncases = ncases+2;

        %% (2) Same state twice, a o a, by name
        f = build_f(Ntot,[1,2],'aa');
        fnew = subs(f,{'a'},ca);
        assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;b}),tol,label+" (2)");
        ncases = ncases+1;

        %% (3) Single-factor monomials, function and functional
        for kind = {'function','functional'}
            f = build_f(Ntot,[1,2],'a',kind{1});
            fnew = subs(f,{'a'},ca);
            assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;b}),tol,label+" (3) "+kind{1});
            ncases = ncases+1;
        end

        % ... also with a Kronecker-type (single-factor) output, for which the
        % evaluation renames the output variables s -> s_1
        f = build_f(Ntot,[2,2],'a');
        fnew = subs(f,{'a'},ca);
        Yref = eval_polyopvar(f,{ca;b});
        Yref.varname = regexprep(Yref.varname,"_1$","");
        assert_close(eval_new(fnew,{}),Yref,tol,label+" (3b)");
        ncases = ncases+1;
        %% (4) Several monomials: those with substituted states change, the rest is kept
        f = build_f(Ntot,[1,2],'a+b');
        fnew = subs(f,{'a'},ca);
        assert_close(eval_new(fnew,{'b',b}),eval_polyopvar(f,{ca;b}),tol,label+" (4)");
        fnew = subs(f,{'a';'b'},[ca;cb]);
        assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (4b)");
        ncases = ncases+2;

        %% (5) Factors by position
        f = build_f(Ntot,[1,2],'ab');
        fnew = subs(f,[1,2],[ca;cb]);
        assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (5)");
        fnew = subs(f,[2,1],[cb;ca]);
        assert_close(eval_new(fnew,{}),eval_polyopvar(f,{ca;cb}),tol,label+" (5b)");
        ncases = ncases+2;

        %% (6) Unsupported cases must error
        f = build_f(Ntot,[1,2],'ab');
        % Partial monomials are supported: the substituted factor remains as a known factor
        % (see also Test_polyopvar_subs_const_partial)
        aa = rand_val(svs,'polynomial');
        for kind = {'pw','kron'}
            if strcmp(kind{1},'pw'),    fp = f;
            else,                       fp = build_f(Ntot,[2,2],'ab');      end
            fnew = subs(fp,{'a'},ca);
            assert_close(eval_new(fnew,{'b',b}),eval_polyopvar(fp,{ca;b}),tol,label+" (6) partial a, "+kind{1});
            fnew = subs(fp,1,ca);
            assert_close(eval_new(fnew,{'b',b}),eval_polyopvar(fp,{ca;b}),tol,label+" (6) partial factor 1, "+kind{1});
            fnew = subs(fp,{'b'},cb);
            assert_close(eval_new(fnew,{'a',aa}),eval_polyopvar(fp,{aa;cb}),tol,label+" (6) partial b, "+kind{1});
            fnew = subs(fp,2,cb);
            assert_close(eval_new(fnew,{'a',aa}),eval_polyopvar(fp,{aa;cb}),tol,label+" (6) partial factor 2, "+kind{1});
            ncases = ncases+4;
        end
        expect_error(@() subs(f,{'a';'b'},ca),'Number of rows');            % wrong number of rows
        expect_error(@() subs(f,{'c'},ca),'states');                        % unknown state
        expect_error(@() subs(f,{'a';'a'},[ca;ca]),'unique');               % repeated state
        expect_error(@() subs(f,{'a';'b'},[ca,ca]),'column');               % not a column
        if Ntot==1
            expect_error(@() subs(f,{'a';'b'},[s2;s2]),'spatial variables');   % wrong variable
        end
        fk = build_f(Ntot,[2,2],'ab');
        expect_error(@() subs(fk,{'a';'b'},[ca;cb]),'tensor product');      % Kronecker output
        ncases = ncases+1;
    end
end

fprintf('Test_polyopvar_subs_const passed (%d substitution cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function f = build_f(N,ot,mono,kind)
% Polynomial in states a,b (both scalar, depending on all N variables):
%   'ab':   f = Kop*(a o b)          'aa':   f = Kop*(a o a)
%   'a':    f = Kop*a                'a+b':  f = Kop1*a + Kop2*b
% as function (pointwise/Kronecker type OT) or, for 'a', functional.
if nargin<4,    kind = 'function';  end
svars = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],N,1));
vars = struct('in',{svars},'out',{svars});
switch mono
    case 'ab'
        f = make_f([1,2],2,N,svars,dom,vars,ot,kind);
    case 'aa'
        f = make_f([1,1],2,N,svars,dom,vars,ot,kind);
    case 'a'
        f = make_f(1,1,N,svars,dom,vars,ot,kind);
    case 'a+b'
        f = make_f(1,1,N,svars,dom,vars,ot,kind) + make_f(2,1,N,svars,dom,vars,ot,kind);
end
end

function f = make_f(Z,d,N,svars,dom,vars,ot,kind)
f = polyopvar();
f.Z = Z;                    f.varname = {'a';'b'};
f.dims.in = [1;1];          f.dims.out = 1;
f.vars.in = svars;          f.dom.in = dom.in;
f.depmat.in = true(2,N);
if strcmp(kind,'function')
    dep = struct('in',true(d,N),'out',true(d,N));
    Kop = rand_tensopvar(ones(d,2),vars,dom,dep,1,ot,2);
    if ~Kop.chkval,     Kop = initialize(Kop);      end
    f.C = tensopmat(Kop);
    f.vars.out = svars;     f.dom.out = dom.out;    f.depmat.out = Kop.deps.out;
else
    Aop = rand_intop([1,1],svars,dom.in,ones(1,N),3);
    f.C = tensopmat(Aop);
    f.vars.out = cell(1,0); f.dom.out = zeros(0,2); f.depmat.out = zeros(1,0);
end
f = initialize(f);
end

function v = rand_val(svs,form)
% Random scalar constant (function) of the given type
if strcmp(form,'double')
    v = randn;
else
    v = polynomial(randn);
    for i=1:numel(svs)
        v = v + randn*svs{i} + randn*svs{i}^2;
    end
end
end

function Y = eval_new(fnew,vals)
% Evaluate fnew, giving values {'name',value,...} to the remaining states
xv = cell(numel(fnew.varname),1);
for k=1:numel(xv)
    xv{k} = polynomial(0);
    for i=1:2:numel(vals)
        if strcmp(fnew.varname{k},vals{i}),     xv{k} = vals{i+1};     end
    end
end
Y = eval_polyopvar(fnew,xv);
end

function assert_deg0(fnew,label)
if any(any(fnew.Z))
    error("Test_polyopvar_subs_const failed (%s): result should be a constant (degree-0) term.",label)
end
end

function assert_close(Ynew,Yref,tol,label)
Yc = polynomial(Yref).coeff;
if isempty(Yc) || max(abs(full(Yc(:))))<1e-6
    error("Test_polyopvar_subs_const (%s): reference value is zero; test would be vacuous.",label)
end
err = polynomial(Ynew-Yref);
if ~isempty(err.coeff) && max(abs(full(err.coeff(:))))>tol
    error("Test_polyopvar_subs_const failed (%s): residual %.3e.",label,max(abs(full(err.coeff(:)))))
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_polyopvar_subs_const failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_polyopvar_subs_const failed: expected an error containing '%s', but none was thrown.",str)
end
