function Test_tensopvar_plus()
% Test_tensopvar_plus.m -- one local function for each case:
%   case_basic              @tensopvar/plus.m. The oracle is
%                           apply_op(A+B,v) against
%                           apply_op(A,v)+apply_op(B,v) on random inputs v
%   case_nested             tensopvar/plus for a nested operator (a sum of
%                           operators with different dependence on the
%                           variables, the operator of one factor of which
%                           is itself a 'tensopvar',
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_tensopvar_plus_nested into this test, as
%   cases;

cases = {@case_basic, @case_nested};
nfail = 0;
for i=1:numel(cases)
    tstart = tic;
    try
        cases{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(cases{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(cases{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_tensopvar_plus: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_tensopvar_plus: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- tests @tensopvar/plus.m. The oracle is
% apply_op(A+B,v) against apply_op(A,v)+apply_op(B,v) on random inputs v.
% otype(1) alternates over {1,2} (pointwise vs Kronecker rows), which
% take different paths in plus.
%

% NOTES
% - otype(2)=2 only: apply_op requires it, so otype(2)==1 cannot be
%   verified through this route.

    ntrials = 20;
    for trial=1:ntrials
    % Set parameters defining tensor-PI operators
    nterms = randi(3);      % number of terms in tensor-PI operator
    d = randi(5);           % degree of the tensor-PI operator
    N = randi(3);           % number of input variables
    M = randi(3);           % number of output variables
    if mod(trial,2)==1
        type = [1,2];       % pointwise rows, Kronecker columns
    else
        type = [2,2];       % Kronecker rows and columns
    end
    dnsty = 0.25;
    % Set parameters defining the function
    p = 2;              % column dimension of function
    pdegv = 2;

    % Matrix dimensions of each factor in the tensor-PI operators
    dims = randi(4,[d,2]);
    if type(1)==1
        m = dims(1,1);
        dims(:,1) = m;
    else
        m = prod(dims(:,1));
    end
    if type(2)==1
        n = dims(1,2);
        dims(:,2) = n;
    else
        n = prod(dims(:,2));
    end
    % Input and output variables
    vars = struct();
    % disjoint 'x'/'y' name pools: input and output variables must not share a name (rand_sopvar rejects one name with two domains)
    vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
    vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
    % Domains
    dom = struct();
    dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
    dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
    % Monomial degrees
    degs = struct();
    degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
    % Dependency arrays
    depmat = struct();
    depmat.in = logical(randi([0,1],d,N));
    depmat.out = logical(randi([0,1],d,M));




    % Generate the random tensor-PI operators and function
    Aop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
    Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);

    v = cell(1,d);
    for j=1:d
        pvars_j = polynomial(Aop.ops{1,j}.vars.in(:));
        n_j = size(Aop.ops{1,j},2);
        v{j} = rand_poly([n_j,p],pvars_j,pdegv,50);
    end

    Cop = Aop + Bop;

    Av = apply_op(Aop,v);
    Bv = apply_op(Bop,v);
    Cv = apply_op(Cop,v);

    tol = 1e-8;
    err_C = cleanpoly(Cv - (Av+Bv),1e-10);
    disp('err_C = ')
    disp(max(max(abs(err_C.C))))
    if max(max(abs(err_C.C)))>tol
        error("Trial %d (type=[%d,%d]): outputs do not match: residual exceeds tolerance.",trial,type(1),type(2))
    end
    end

    fprintf('Test_tensopvar_plus passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
end

%% --------------------------------------------------------------------------
function case_nested()
% case_nested -- tests tensopvar/plus for a nested
% operator (a sum of operators with different dependence on the
% variables, the operator of one factor of which is itself a 'tensopvar',
% d=1) plus an operator with several factors (d>1), in either order and
% in any bracketing of a sum of three or four operators.
%
% Verification is semantic and independent of the routine under test: the
% operators are applied to the factors x, 1, y of a monomial x*m*y
% (apply_op, which takes one input for each factor, counting those of
% nested operators) and compared with the polynomial
% x*(m_a+m_b+...)*y.
% Parts:
%   (1) sums of three operators in all orders and bracketings;
%   (2) sums of four operators, left-to-right, and the sum of a nested
%       operator with a nested operator of the same structure;
%   (3) sums that do not involve a nested operator (single factors, equal
%       structure) are as for non-nested operators.
%

    tol = 1e-8;
    ncases = 0;
    pvar s1 s2
    m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;      m4 = 1+s1+s2^2;
    polyopvar x(s1,s2) y(s1,s2)
    Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1});
    xv = 1+s1;          yv = 2+s2-s1*s2;
    ms = {m1,m3,m4,m2};                 % multipliers, depending on s1, s2, s1 and s2, s1 and s2
    T = cell(1,4);
    for k=1:4
        T{k} = build_op({'x',ms{k},'y'},Id);
    end
    leaves = {xv,1,yv};

    %% (1) Three operators, in all orders and bracketings
    orders = perms(1:3);
    for r=1:size(orders,1)
        o = orders(r,:);
        Yref = xv*(ms{o(1)}+ms{o(2)}+ms{o(3)})*yv;
        left = (T{o(1)}+T{o(2)})+T{o(3)};
        right = T{o(1)}+(T{o(2)}+T{o(3)});
        check_nested(left,sprintf('(1.%d) left',r));
        assert_close(apply_op(left,leaves),Yref,tol,sprintf('(1.%d) left',r));
        assert_close(apply_op(right,leaves),Yref,tol,sprintf('(1.%d) right',r));
        ncases = ncases+2;
    end

    %% (2) Four operators
    Yref = xv*(m1+m3+m4+m2)*yv;
    left = ((T{1}+T{2})+T{3})+T{4};
    assert_close(apply_op(left,leaves),Yref,tol,'(2.1) left-to-right');
    mixed = (T{1}+(T{2}+T{3}))+T{4};
    assert_close(apply_op(mixed,leaves),Yref,tol,'(2.2) mixed brackets');
    % nested + nested: (Ta+Tb) + (Tc+Td), and the sum of a nested operator with itself
    bal = (T{1}+T{2})+(T{3}+T{4});
    assert_close(apply_op(bal,leaves),Yref,tol,'(2.3) balanced');
    Tab = T{1}+T{2};
    assert_close(apply_op(Tab+Tab,leaves),2*xv*(m1+m3)*yv,tol,'(2.4) nested+itself');
    ncases = ncases+4;

    %% (3) Without a nested operator: unchanged
    % equal structure and dependence: factors concatenated as terms, value = sum
    S = T{1}+T{1};
    assert_close(apply_op(S,leaves),2*xv*m1*yv,tol,'(3.1) equal structure');
    % single factors: the factors are added
    vars = struct('in',{{'s1'}},'out',{{'s1'}});
    dom = struct('in',[0,1],'out',[0,1]);
    M1 = mat2tensopvar(1+s1,vars,dom);          M2 = mat2tensopvar(2-s1^2,vars,dom);
    v = 1+2*s1;
    Y = apply_op(M1+M2,v);
    Yref = (1+s1)*v+(2-s1^2)*v;
    Y.varname = regexprep(Y.varname,'_1$','');
    assert_close(Y,Yref,tol,'(3.2) single factors');
    % a single factor plus an operator with two factors, of equal dependence: the sum cannot act on a
    % common input, but the second factor of B must not be dropped silently
    A1 = mat2tensopvar(1+0*s1,vars,dom);
    mR = mat2tensopvar(2+s1,struct('in',{cell(1,0)},'out',{{'s1'}}),struct('in',zeros(0,2),'out',[0,1]));
    B2 = otimes(A1,mR,[1,2]);
    C = A1+B2;
    if ~isequal(size(C.ops),[2,1]) || ~isa(C.ops{2,1},'tensopvar') || num_factors(C.ops{2,1})~=2
        error("case_nested failed (3.3): the factors of the second operator were dropped.")
    end
    ncases = ncases+3;

    fprintf('case_nested passed (%d cases).\n',ncases)
end

%% ---- helpers --------------------------------------------------------------

function K = build_op(facs,Id)
% Tensor-PI operator acting on the monomial: the identity operators of the states for 'x','y', multipliers
% from R for polynomials, in a pointwise product (otype [1,2])
K = [];
for j=1:numel(facs)
    if ischar(facs{j})
        Fj = Id.(facs{j});
    else
        vars = cell(1,0);
        if ismember('s1',facs{j}.varname),  vars = [vars,{'s1'}];   end
        if ismember('s2',facs{j}.varname),  vars = [vars,{'s2'}];   end
        Fj = mat2tensopvar(facs{j},struct('in',{cell(1,0)},'out',{vars}),struct('in',zeros(0,2),'out',repmat([0,1],numel(vars),1)));
    end
    if isempty(K),  K = Fj;
    else,           K = otimes(K,Fj,[1,2]);     end
end
end

function check_nested(Top,label)
if ~any(cellfun('isclass',Top.ops(:),'tensopvar'))
    error("case_nested failed %s: setup: the sum should be a nested operator.",label)
end
end

function assert_close(Y,Yref,tol,label)
% Output variables of Kronecker products are named s_1,..., renamed back
Y = polynomial(Y);      Yref = polynomial(Yref);
Y.varname = regexprep(Y.varname,'_1$','');
D = Y-Yref;
if ~isequal(size(D),size(Yref))
    error("case_nested failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol
    error("case_nested failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient))
    error("case_nested %s: reference value is zero; test would be vacuous.",label)
end
end
