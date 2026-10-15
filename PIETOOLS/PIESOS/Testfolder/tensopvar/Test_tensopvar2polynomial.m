% Test_tensopvar2polynomial.m -- tests tensopvar2polynomial, converting a
% 'tensopvar' that acts on a constant monomial (no input variables, all
% factors multipliers from R) to the 'polynomial' P with
%   Top(x1,...,xd) = P*(x1 o ... o xd),
% o the Kronecker product of the inputs. The expected values are the
% explicit matrices of the multipliers, and the action of the operator on
% random inputs, apply_op, which does not use the conversion.
% Parts:
%   (1) a single factor: scalar, vector and matrix multipliers,
%       constants, a 'sopvar';
%   (2) several factors and terms with pointwise rows (random), compared
%       with apply_op;
%   (3) a Kronecker product of the rows, for factors in different
%       variables, compared with apply_op;
%   (4) nested operators, and pointwise columns;
%   (5) operators that cannot be converted: IS_OK is false and P empty, or
%       an error without the output IS_OK.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-9;
pvar s1 s2 s3
ncases = 0;
dom1 = [0,1];
mk = @(M,vars) mat2tensopvar(M,struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat(dom1,numel(vars),1)));

%% (1) A single factor
Ms = {  1+s1^2,             {'s1'}
        [1+s1;2],           {'s1'}
        [s1,s2;3,1-s1*s2],  {'s1','s2'}
        polynomial(3),      cell(1,0)
        [2;5],              cell(1,0)           };
for k=1:size(Ms,1)
    M = Ms{k,1};
    [P,is_ok] = tensopvar2polynomial(mk(M,Ms{k,2}));
    if ~is_ok || ~isa(P,'polynomial')
        error("Test_tensopvar2polynomial failed (1.%d): not converted.",k)
    end
    assert_close(P,M,tol,sprintf('(1.%d)',k));
    ncases = ncases+1;
end
P = tensopvar2polynomial(mk([1+s1;2],{'s1'}));          % one output: error if not converted, so no error here
assert_close(P,[1+s1;2],tol,'(1) one output');
% a 'sopvar'
Ps = tensopvar2polynomial(mat2sopvar([s1;1],struct('in',{cell(1,0)},'out',{{'s1'}}),struct('in',zeros(0,2),'out',dom1)));
assert_close(Ps,[s1;1],tol,'(1) sopvar');
ncases = ncases+2;

%% (2) Several factors and terms, pointwise rows, Kronecker columns, compared with apply_op
for rep=1:24
    rng(rep);
    d = 1+mod(rep,3);               nterms = 1+mod(floor(rep/3),3);         M = 1+mod(rep,2);
    nr = 1+randi(2);                nc = randi(3,d,1);
    dims = [nr*ones(d,1),nc];
    vars = struct('in',{cell(1,0)},'out',{arrayfun(@(i) ['s',num2str(i)],1:M,'UniformOutput',false)});
    dom = struct('in',zeros(0,2),'out',repmat(dom1,M,1));
    depmat = struct('in',false(d,0),'out',rand(d,M)>0.3);
    Top = rand_tensopvar(dims,vars,dom,depmat,nterms,[1,2],2,[]);
    [P,is_ok] = tensopvar2polynomial(Top);
    if ~is_ok || ~isequal(size(P),[nr,prod(nc)])
        error("Test_tensopvar2polynomial failed (2.%d): not converted, or wrong size.",rep)
    end
    xs = arrayfun(@(j) randn(nc(j),1),1:d,'UniformOutput',false);
    xx = 1;     for j=1:d,  xx = kron(xx,xs{j});    end
    assert_close(apply_op(Top,xs),P*xx,tol,sprintf('(2.%d)',rep));
    ncases = ncases+1;
end

%% (3) Kronecker product of the rows: factors in different variables
for rep=1:6
    rng(50+rep);
    d = 2+mod(rep,2);               nterms = 1+mod(rep,2);
    nr = randi(2,d,1);              nc = randi(2,d,1);
    dims = [nr,nc];
    vars = struct('in',{cell(1,0)},'out',{arrayfun(@(i) ['s',num2str(i)],1:d,'UniformOutput',false)});
    dom = struct('in',zeros(0,2),'out',repmat(dom1,d,1));
    depmat = struct('in',false(d,0),'out',eye(d)>0);        % factor j depends on variable sj only
    Top = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],2,[]);
    [P,is_ok] = tensopvar2polynomial(Top);
    if ~is_ok || ~isequal(size(P),[prod(nr),prod(nc)])
        error("Test_tensopvar2polynomial failed (3.%d): not converted, or wrong size.",rep)
    end
    xs = arrayfun(@(j) randn(nc(j),1),1:d,'UniformOutput',false);
    xx = 1;     for j=1:d,  xx = kron(xx,xs{j});    end
    Y = apply_op(Top,xs);       Y.varname = regexprep(Y.varname,'_1$','');     % apply_op names the variables s1_1
    assert_close(Y,P*xx,tol,sprintf('(3.%d)',rep));
    ncases = ncases+1;
end

%% (4) Nested operators, and pointwise columns
Ta = otimes(mk([1+s1;2],{'s1'}),mk([s2,3;1,s2],{'s2'}),[1,2]);          % depends on s1, s2 in different factors
Tb = otimes(mk([s1;s1^2],{'s1','s2'}),mk([2,s1;s1,1],{'s1','s2'}),[1,2]);         % both factors on s1, s2: the sum is nested
Tab = Ta+Tb;
if ~isa(Tab.ops{1,1},'tensopvar')
    error("Test_tensopvar2polynomial failed (4.1): setup: the sum should be a nested operator.")
end
[P,is_ok] = tensopvar2polynomial(Tab);
if ~is_ok
    error("Test_tensopvar2polynomial failed (4.1): nested operator not converted.")
end
Pa = tensopvar2polynomial(Ta);      Pb = tensopvar2polynomial(Tb);
assert_close(P,Pa+Pb,tol,'(4.1)');
x1 = 2;     x2 = [1;-2];
Y = apply_op(Tab,{{x1,x2}});         Y.varname = regexprep(Y.varname,'_1$','');
assert_close(Y,P*kron(x1,x2),tol,'(4.1) values');
% pointwise columns (otype(2)=1): P(r,c) = F1(r,c)*F2(r,c)
Tp = otimes(mk([1+s1,2;s1,3],{'s1'}),mk([s2,1;2,s2],{'s2'}),[1,1]);
P = tensopvar2polynomial(Tp);
assert_close(P,[1+s1,2;s1,3].*[s2,1;2,s2],tol,'(4.2)');
ncases = ncases+3;

%% (5) Operators that cannot be converted
vars_in = struct('in',{{'s1'}},'out',{{'s1'}});         dom_in = struct('in',dom1,'out',dom1);
Tin = mat2tensopvar(1+s1,vars_in,dom_in);               % multiplier on L2: acts on a function space
Tfn = mat2tensopvar(2,struct('in',{{'s1'}},'out',{cell(1,0)}),struct('in',dom1,'out',zeros(0,2)));   % a functional
Tsh = otimes(mk([1+s1;2],{'s1'}),mk([s1;1],{'s1'}),[2,2]);   % Kronecker product in the same variable
Kint = sopvar2tensopvar({intop({1},{[0,2]},{'s1'},[0,1],true)});    % an 'intop' as factor
cases5 = {Tin,'function space'; Tfn,'function space'; Tsh,'cannot be converted'; Kint,'function space'; 5,'tensopvar'};
for k=1:size(cases5,1)
    [P,is_ok] = tensopvar2polynomial(cases5{k,1});
    if is_ok || ~isempty(P)
        error("Test_tensopvar2polynomial failed (5.%d): an operator that cannot be converted was converted.",k)
    end
    try
        tensopvar2polynomial(cases5{k,1});
        error("Test_tensopvar2polynomial failed (5.%d): no error without the output IS_OK.",k)
    catch ME
        if ~contains(ME.message,cases5{k,2})
            error("Test_tensopvar2polynomial failed (5.%d): unexpected error message '%s' (expected to contain '%s').",k,ME.message,cases5{k,2})
        end
    end
    ncases = ncases+1;
end

fprintf('Test_tensopvar2polynomial passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function assert_close(Y,Yref,tol,label)
Yref = polynomial(Yref);
D = polynomial(Y)-Yref;
if ~isequal(size(D),size(Yref))
    error("Test_tensopvar2polynomial failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol
    error("Test_tensopvar2polynomial failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient))
    error("Test_tensopvar2polynomial %s: reference value is zero; test would be vacuous.",label)
end
end
