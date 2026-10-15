% Test_tensopvar_apply_op_nested.m -- tests apply_op of a 'tensopvar' with
% nested factors (a 'tensopvar' as factor of another), applied to a cell
% with one input for each factor, counting those of the nested operators
% (num_factors(Cop)), as apply_poly/subs pass them for the operator acting
% on a monomial of a 'polyopvar', as well as to the 1 x d cell with the
% inputs of nested factors as cells.
%
% Such operators arise from the sum of monomials with the same Z and the
% same number of factors, with factors that depend on different variables
% (tensopvar/plus), and as the pointwise product of operators in a
% Kronecker product (nonlinear terms of a PIE).
% Verification is semantic and independent of the routine under test:
%   (1) 'polyopvar': the sum of monomials x*m1*y + x*m2*y with the m_k in
%       different variables, evaluated by apply_poly, vs. the product of
%       the polynomials x*(m1+m2)*y;
%   (2) 'tensopvar': the output of apply_op for the flat cell of leaf
%       inputs vs. (a) the value for the 1 x d form of the inputs, with
%       the inputs of nested factors as cells (built here by
%       nest_leaves), and (b) apply_nested, a recursion over the nested
%       operators that takes the inputs of the leaf operators one by one
%       from the flat list;
%   (3) operators of which the nesting differs between the terms
%       (apply_nested), which the 1 x d form cannot represent;
%   (4) wrong numbers or dimensions of inputs are rejected (the dimension
%       of the input of a factor after a nested one is checked against the
%       right factor), and non-nested operators are unaffected.
%
% Sums of three or more monomials are built from the left and from the
% right (cases 1.x).
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Kronecker products of operators nested in Kronecker products (cases
%   2.1, 2.2) are only compared to the 1 x d form: apply_op names the
%   output variables one per factor of the outer operator (s_1,s_2), where
%   apply_nested numbers each leaf operator (s_1,s_2,s_3);
% - Substitution of some, not all, states in a nested operator: see
%   Test_tensopvar_mtimes_nested.

tol = 1e-8;
ncases = 0;
pvar s s1 s2

%% (1) polyopvar: nested sums of monomials with known factors, evaluated by apply_poly
m1 = 1-s1^2;        m2 = 2+s2+s1*s2;        m3 = 3-s2;      m4 = 1+s1+s2^2;
polyopvar x(s1,s2) y(s1,s2)
Id = struct('x',x.C.ops{1,1},'y',y.C.ops{1,1});
xv = 1+s1;          yv = 2+s2-s1*s2;
% Each row: Z of the monomials, factors, coefficients, for which the sum is nested
specs = {   [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'}},                  [2,-1]
            [0,1,2],    {{m1,'x','y'},{m3,'x','y'}},                  [1,1]
            [1,2,0],    {{'x','y',m1},{'x','y',m3}},                  [1,1]
            [1,0],      {{'x',m1},{'x',m3}},                          [1,1]
            [1,0,0,2],  {{'x',m1,m3,'y'},{'x',m3,m1,'y'}},            [1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'},{'x',m4,'y'}},     [1,1,1]
            [1,0,2],    {{'x',m1,'y'},{'x',m3,'y'},{'x',m4,'y'}},     [1,-2,3]   };
for k=1:size(specs,1)
    Zrow = specs{k,1};      facs = specs{k,2};      cf = specs{k,3};
    % Sum from the left, (f1+f2)+f3, and from the right, f1+(f2+f3)
    Yref = 0;       fl = cell(1,numel(facs));
    for l=1:numel(facs)
        fl{l} = cf(l)*build_f(Zrow,facs{l},Id);
        Yref = Yref + cf(l)*prod_facs(Zrow,facs{l},xv,yv);
    end
    f = fl{end};
    for l=numel(fl)-1:-1:1
        f = fl{l} + f;
    end
    fleft = fl{1};
    for l=2:numel(fl)
        fleft = fleft + fl{l};
    end
    for f = {combine_terms(f),combine_terms(fleft)}
        f = f{1};
        Top = f.C.ops{1,1};
        if ~isa(Top,'tensopvar') || ~any(cellfun('isclass',Top.ops(:),'tensopvar'))
            error("Test_tensopvar_apply_op_nested failed (1.%d): setup: the operator should be nested.",k)
        end
        assert_close(apply_poly(f,{xv,yv}),Yref,tol,sprintf('(1.%d)',k));
        ncases = ncases+1;
    end
end

%% (2) tensopvar: flat inputs vs the 1 x d form of the inputs, and vs apply_nested
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
dep1 = struct('in',true,'out',true);
pd = struct('in',1,'out',1);
mult = @() mat2tensopvar(randn+randn*s+randn*s^2,vars,dom).ops{1,1};      % multiplier on L2[0,1]
intl = @() rand_tensopvar([1,1],vars,dom,dep1,1,[2,2],pd,0.5).ops{1,1};    % random 3-PI operator
for variant=1:6
    rng(100+variant);
    switch variant
    case 1      % {{N1 o S1} o N2}, N1 = P1*P2 and N2 = P3*P4 pointwise products, S1 an integral: 5 leaves
        N1 = pointwise(mult(),mult());     S1 = intl();      N2 = pointwise(mult(),mult());
        R12 = kron_op({N1,S1});            Bop = kron_op({R12,N2});
    case 2      % {N1 o {S1 o N2}}, with an integral in N2: 5 leaves
        N1 = pointwise(mult(),mult());     S1 = mult();      N2 = pointwise(mult(),intl());
        R12 = kron_op({S1,N2});            Bop = kron_op({N1,R12});
    case 3      % one nested pointwise product per term
        Bop = kron_op({pointwise(mult(),mult());pointwise(mult(),mult())});
    case 4      % pointwise product of three operators as the first factor, then an integral, two terms
        Bop = kron_op({pointwise(pointwise(mult(),mult()),mult()),intl(); ...
                       pointwise(pointwise(mult(),mult()),mult()),intl()});
    case 5      % a nested operator with a single factor: as many inputs as factors
        Bop = kron_op({kron_op({mult()})});
    case 6      % nested operator with several terms, each a pointwise product, two factors in the outer operator
        N1 = pointwise_terms(mult(),mult(),mult(),mult());
        Bop = kron_op({N1,mult()});
    end
    nl = num_factors(Bop);
    pvars = polynomial(Bop.vars.in(:));
    leaves = cell(1,nl);
    for l=1:nl
        leaves{l} = rand_poly([1,1],pvars,2,10);
    end
    Yflat = apply_op(Bop,leaves);
    Ystruct = apply_op(Bop,nest_leaves(Bop,leaves));
    assert_close(Yflat,Ystruct,tol,sprintf('(2.%d) vs 1 x d form',variant));
    if variant>=3
        Yref = apply_nested(Bop,leaves,1,zeros(1,numel(Bop.vars.out)),Bop.vars.out);
        assert_close(Yflat,Yref,tol,sprintf('(2.%d) vs apply_nested',variant));
    end
    ncases = ncases+1;
end

%% (3) The nesting differs between the terms: {N(2 leaves), S} and {S, N(2 leaves)}
% multipliers from R, so that the dependency of the factors on the variables is the same in each term
mR = @() mat2tensopvar(randn+randn*s+randn*s^2,struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1])).ops{1,1};
for variant=1:4
    rng(200+variant);
    switch variant
    case 1
        Bop = kron_op({pointwise(mR(),mR()),mR(); mR(),pointwise(mR(),mR())});
    case 2
        Bop = kron_op({pointwise(mR(),mR()),mR(); mR(),pointwise(mR(),mR()); pointwise(mR(),mR()),mR()});
    case 3
        Bop = kron_op({pointwise(mR(),mR()),mR(),mR(); mR(),mR(),pointwise(mR(),mR())});
    case 4  % nesting of depth two, in different factors
        Bop = kron_op({pointwise(pointwise(mR(),mR()),mR()),mR(); mR(),pointwise(pointwise(mR(),mR()),mR())});
    end
    nl = num_factors(Bop);
    leaves = num2cell(ones(1,nl));
    leaves = cellfun(@(c) c*(1+rand),leaves,'UniformOutput',false);      % column dimension 1: any scalar
    Yflat = apply_op(Bop,leaves);
    Yref = apply_nested(Bop,leaves,1,zeros(1,numel(Bop.vars.out)),Bop.vars.out);
    assert_close(Yflat,Yref,tol,sprintf('(3.%d)',variant));
    ncases = ncases+1;
end

%% (4) Wrong numbers or dimensions of inputs, and operators without nesting
% a factor after the nested one with a column dimension of 2: its input is the third of the flat list
mR2 = mat2tensopvar([1+s,2*s^2],struct('in',{cell(1,0)},'out',{{'s'}}),struct('in',zeros(0,2),'out',[0,1])).ops{1,1};
Bop = kron_op({pointwise(mR(),mR()),mR2});
leaves = {1.3,0.7,[1;-2]};
assert_close(apply_op(Bop,leaves),apply_nested(Bop,leaves,1,zeros(1,1),Bop.vars.out),tol,'(4) column dimension 2');
expect_error(@() apply_op(Bop,{1.3,0.7,1}),'Column dimension');
expect_error(@() apply_op(Bop,{1.3,[1;-2],0.7}),'Column dimension');
Bop = kron_op({pointwise(mR(),mR()),mR()});         % 2 factors, 3 leaves
expect_error(@() apply_op(Bop,{1,1}),'1xd cell');       % d inputs: the 1 x d form, and a number is not the cell of a nested factor
expect_error(@() apply_op(Bop,{1,1,1,1}),'input for each factor');
expect_error(@() apply_op(Bop,{1}),'input for each factor');
% the terms have different numbers of leaves: {N(2), S} and {N(3), S}
Bad = kron_op({pointwise(mR(),mR()),mR(); pointwise(pointwise(mR(),mR()),mR()),mR()});
expect_error(@() apply_op(Bad,{1,1,1}),'number of inputs must match');
% a 1 x d cell for an operator without nesting, and a bare input for a single factor, as before
T = kron_op({mult(),mult()});
v1 = rand_poly([1,1],polynomial({'s'}),2,10);       v2 = rand_poly([1,1],polynomial({'s'}),2,10);
Y = apply_op(T,{v1,v2});
Ya = apply_sopvar(T.ops{1,1},v1);       Yb = apply_sopvar(T.ops{1,2},v2);
Ya.varname(strcmp(Ya.varname,'s')) = {'s_1'};           Yb.varname(strcmp(Yb.varname,'s')) = {'s_2'};
assert_close(Y,Ya*Yb,tol,'(4) no nesting');
T1 = kron_op({mult()});
Y1 = apply_op(T1,v1);       Y1c = apply_op(T1,{v1});
assert_close(Y1,Y1c,tol,'(4) single factor');
expect_error(@() apply_op(T,v1),'1xd cell');
expect_error(@() apply_op(T,{v1,v2,v1}),'1xd cell');
ncases = ncases+11;

fprintf('Test_tensopvar_apply_op_nested passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function N = pointwise(A,B)
% Pointwise product of the outputs of the operators A and B
N = tensopvar();        N.otype = [1,2];        N.ops = {A,B};      N = initialize(N);
end

function N = kron_op(ops)
% Kronecker product of the outputs of the operators in each row of the cell ops, summed over the rows
N = tensopvar();        N.otype = [2,2];        N.ops = ops;        N = initialize(N);
end

function N = pointwise_terms(A1,B1,A2,B2)
% Sum of two pointwise products, a nested operator with two terms
N = tensopvar();        N.otype = [1,2];        N.ops = {A1,B1; A2,B2};     N = initialize(N);
end

function X = nest_leaves(R,leaves)
% The 1 x d cell of inputs to the first term of R, with the inputs of nested factors as cells
pos = 1;        d = size(R.ops,2);      X = cell(1,d);
for j=1:d
    F = R.ops{1,j};
    if isa(F,'tensopvar')
        n = num_factors(F);
        X{j} = nest_leaves(F,leaves(pos:pos+n-1));
        pos = pos+n;
    else
        X{j} = leaves{pos};
        pos = pos+1;
    end
end
end

function [Y,pos,nslots] = apply_nested(R,leaves,pos,off,vars_out)
% Apply the (nested) tensor-PI operator R to the states in leaves(pos:end), one for each of the 'sopvar'
% operators in a term of R, in order, each term acting on the same states. The output variable vars_out{k}
% is replaced by vars_out{k}_(off(k)+1), vars_out{k}_(off(k)+2),... for the leaf operators in a Kronecker
% product of the outputs, whereas the outputs of the factors in a pointwise product share the variable.
M = numel(vars_out);
if isa(R,'sopvar')
    Y = apply_sopvar(R,leaves{pos});
    pos = pos+1;
    nslots = zeros(1,M);
    for k=1:M
        if ismember(vars_out{k},R.vars.out)
            is_v = strcmp(Y.varname,vars_out{k});
            Y.varname(is_v) = {[vars_out{k},'_',num2str(off(k)+1)]};
            nslots(k) = 1;
        end
    end
    return
end
Y = 0;
pos0 = pos;
for i=1:size(R.ops,1)
    Yi = 1;     pos = pos0;     off_j = off;    nslots = zeros(1,M);
    for j=1:size(R.ops,2)
        [Yij,pos,nslots_j] = apply_nested(R.ops{i,j},leaves,pos,off_j,vars_out);
        if R.otype(1)==1
            Yi = Yi.*Yij;
            nslots = max(nslots,nslots_j);
        else
            Yi = kron(Yi,Yij);
            off_j = off_j+nslots_j;
            nslots = nslots+nslots_j;
        end
    end
    Y = Y+Yi;
end
end

function f = build_f(Zrow,facs,Id)
% polyopvar with monomial Z=Zrow in the states x,y, acted on by the product of the factors facs: the
% identity operators of the states for 'x','y', multipliers from R for polynomials
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
f = polyopvar();
f.varname = {'x';'y'};      f.Z = sparse(Zrow);     f.C = tensopmat(K);
f.dims.in = [1;1];
f.vars.in = {'s1','s2'};    f.dom.in = [0,1;0,1];   f.depmat.in = true(2,2);
f = initialize(f);
end

function Y = prod_facs(Zrow,facs,xv,yv)
% Product of the values of the factors, as polynomials: x -> xv, y -> yv, polynomial -> itself
Y = 1;
for j=1:numel(facs)
    if ischar(facs{j})
        if strcmp(facs{j},'x'),     Y = Y*xv;
        else,                       Y = Y*yv;       end
    else
        Y = Y*facs{j};
    end
end
end

function assert_close(Y,Yref,tol,label)
% Output variables of Kronecker products are named s_1,s_2,..., and compared as they are named
Y = polynomial(Y);      Yref = polynomial(Yref);
D = Y-Yref;
if ~isequal(size(D),size(Yref))
    error("Test_tensopvar_apply_op_nested failed %s: size %s, expected %s.",label,mat2str(size(Y)),mat2str(size(Yref)))
end
err = full(max([0;abs(nonzeros(D.coefficient))]));
if err>tol
    error("Test_tensopvar_apply_op_nested failed %s: residual %.3e.",label,err)
end
if isempty(nonzeros(Yref.coefficient))
    error("Test_tensopvar_apply_op_nested %s: reference value is zero; test would be vacuous.",label)
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_tensopvar_apply_op_nested failed: unexpected error message '%s' (expected to contain '%s').",ME.message,str)
    end
    return
end
error("Test_tensopvar_apply_op_nested failed: expected an error containing '%s', but none was thrown.",str)
end
