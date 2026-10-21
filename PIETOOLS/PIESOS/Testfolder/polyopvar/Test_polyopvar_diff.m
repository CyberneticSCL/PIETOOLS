% TEST_POLYOPVAR_DIFF
% Tests @polyopvar/diff, the derivative of a distributed polynomial P(x)
% with respect to its OUTPUT variables. P has states a,b (vector-valued)
% and monomials a and a o b, with pointwise outputs from 1-3 variables, and
% is also tested for two rows (vertcat). Oracle: the test-only evaluator
% eval_polyopvar (which applies the operators of P with apply_op, and
% shares no code with 'diff') evaluated at polynomial states, against the
% derivative of the evaluation, taken by 'polynomial' calculus, stacked
% over the rows of d.
%
% Run pietools_path_update first.

% DJ (CC), 10/20/2026: Initial coding;

rng(31);
tol = 1e-8;
npass = 0;      nfail = 0;      n_nontriv = 0;      n_err = 0;
dom = [0,1];

for trial=1:20
    N = randi(3);                                   % number of spatial variables
    vn = arrayfun(@(i) ['s',num2str(i)],1:N,'UniformOutput',false);
    sv = polynomial(vn);
    ma = randi(2);      mb = randi(2);
    nrows = randi(2);
    f = build_f(vn,ma,mb);
    if nrows==2
        f = vertcat(f,build_f(vn,ma,mb));
    end
    xa = rand_poly([ma,1],sv,2,5);      xb = rand_poly([mb,1],sv,2,5);
    Y = eval_polyopvar(f,{xa;xb});
    % Differentiate with respect to the first variable, the last, both, rows stacked
    s_first = polynomial(vn(1));        s_last = polynomial(vn(N));
    cases = {{vn(1),1,       @(Y) diff(Y,s_first)}, ...
             {vn(N),1,       @(Y) diff(Y,s_last)}, ...
             {vn([N,1]),[1,1], @(Y) diff(diff(Y,s_first),s_last)}, ...
             {vn([N,1]),[1,0;0,1], @(Y) [diff(Y,s_last);diff(Y,s_first)]}, ...
             {vn(1),2,       @(Y) diff(diff(Y,s_first),s_first)}};
    if N==1
        cases = cases([1,2,4,5]);                   % [N,1] would repeat the variable
    end
    for c_i=1:numel(cases)
        vc = cases{c_i}{1};     dc = cases{c_i}{2};     ref_fun = cases{c_i}{3};
        is_second = isequal(dc,2);
        try
            df = diff(f,vc,dc);
        catch ME
            if is_second && contains(ME.message,"cannot be expressed as a PI operator")
                n_err = n_err+1;
                continue
            else
                rethrow(ME)
            end
        end
        Yref = ref_fun(Y);
        dY = eval_polyopvar(df,{xa;xb});
        n_nontriv = n_nontriv + (norm_poly(Yref)>tol);
        [npass,nfail] = tally(norm_poly(dY-Yref)<tol && isequal(size(dY),size(Yref)) && df.chkval ...
            && isequal(size(df),[size(Yref,1),1]) && isequal(df.Z,f.Z),npass,nfail, ...
            sprintf('trial %d (N=%d,rows=%d), case %d',trial,N,nrows,c_i));
    end
end

%% Edge cases
vn = {'s1','s2'};       sv = polynomial(vn);        s1 = polynomial({'s1'});
f = build_f(vn,2,1);
xa = rand_poly([2,1],sv,2,5);      xb = rand_poly([1,1],sv,2,5);
Y = eval_polyopvar(f,{xa;xb});
% zero order: the polynomial itself, stacked
df = diff(f,{'s1'},[0;0]);
[npass,nfail] = tally(norm_poly(eval_polyopvar(df,{xa;xb})-[Y;Y])<tol,npass,nfail,'zero order, stacked');
% variable that is not an output: zero
df = diff(f,{'s9'},1);
[npass,nfail] = tally(norm_poly(eval_polyopvar(df,{xa;xb}))<tol && isequal(size(df),size(f)),npass,nfail,'non-output variable');
df = diff(f,{'s9','s1'},[1,0;0,1]);
[npass,nfail] = tally(norm_poly(eval_polyopvar(df,{xa;xb})-[0*Y;diff(Y,s1)])<tol,npass,nfail,'zero and nonzero rows');
% state variables are not supported
[npass,nfail] = tally(throws(@() diff(f,{'a'},1),"state variables"),npass,nfail,'state variable error');
[npass,nfail] = tally(throws(@() diff(f,{'s1'},-1),"nonnegative"),npass,nfail,'negative order');
% multiplier in the output variable: not a PI derivative, so an error
fm = build_f(vn,1,1,true);
[npass,nfail] = tally(throws(@() diff(fm,{'s1'},1),"cannot be expressed as a PI operator"),npass,nfail,'multiplier error');

fprintf('\nTest_polyopvar_diff: %d passed, %d failed (%d nonzero comparisons, %d expected second-order errors).\n',npass,nfail,n_nontriv,n_err);
if nfail>0
    error('Test_polyopvar_diff failed.');
end


%% Helpers
function f = build_f(vn,ma,mb,add_mult)
% f = C*Z(x), Z = [1,0; 1,2], i.e. monomials a and a o b, states a (size
% ma) and b (size mb) on all of the variables vn, scalar-row-block output
% of dimension 2 (pointwise product of the factors), sopvar factors without
% multiplier parts in the variables, unless add_mult.
if nargin<4,    add_mult = false;   end
N = numel(vn);
vars = struct('in',{vn},'out',{vn});
dep1 = struct('in',true(1,N),'out',true(1,N));
dep2 = struct('in',true(2,N),'out',true(2,N));
T1 = rand_tensopvar([2,ma],vars,[0,1],dep1,2,[1,2],2,0.3);
T2 = rand_tensopvar([2,ma;2,mb],vars,[0,1],dep2,2,[1,2],2,0.3);
T1 = zero_mult_all(T1,vn);      T2 = zero_mult_all(T2,vn);
if add_mult
    T1 = add_mult_all(T1);
end
f = polyopvar();
f.Z = [1,0;1,2];            f.varname = {'a';'b'};
f.dims.in = [ma;mb];
f.vars.in = vn;             f.dom.in = repmat([0,1],N,1);
f.depmat.in = true(2,N);
f.C = tensopmat({T1,T2});
f = initialize(f);
end

function e = norm_poly(p)
if isa(p,'double')
    e = max([0;abs(p(:))]);
    return
end
p = cleanpoly(polynomial(p),1e-12);
if isempty(p.C) || nnz(p.C)==0
    e = 0;
else
    e = max(abs(p.C(:)));
end
end

function [npass,nfail] = tally(ok,npass,nfail,name)
if ok
    npass = npass+1;
else
    nfail = nfail+1;    fprintf('FAIL: %s\n',name);
end
end

function tf = throws(f,msg)
tf = false;
try
    f();
catch ME
    tf = contains(ME.message,msg);
end
end

function Cop = zero_mult_all(Cop,vn)
% Zero all parameters of all factors with a multiplier in one of the variables vn
for l=1:numel(Cop.ops)
    P = Cop.ops{l};
    vS3 = intersect(P.vars.in,P.vars.out);
    params = P.params;
    for k=1:numel(params)
        g = gamma_of_cell(k,numel(vS3));
        if any(g(ismember(vS3,vn))==1)
            params{k} = sparse(size(params{k},1),size(params{k},2));
        end
    end
    Cop.ops{l} = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
end
Cop = initialize(Cop);
end

function Cop = add_mult_all(Cop)
% Put a nonzero constant in every multiplier parameter of the first factor
P = Cop.ops{1};
vS3 = intersect(P.vars.in,P.vars.out);
params = P.params;
for k=1:numel(params)
    g = gamma_of_cell(k,numel(vS3));
    if all(g==1)
        params{k}(1,1) = 1;
    end
end
Cop.ops{1} = sopvar(params,P.vars,P.ZL,P.ZR,P.dom,P.dims);
Cop = initialize(Cop);
end
