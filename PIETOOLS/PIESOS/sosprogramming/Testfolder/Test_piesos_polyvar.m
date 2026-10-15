function Test_piesos_polyvar()
% TEST_PIESOS_POLYVAR tests, using assertions, piesos_polyvar,
%   [prog,p] = piesos_polyvar(prog,Zx,opts),
% which declares p = Cop*Zx, with Zx a basis of nZ distributed monomials
% (output rows), and Cop = [Kop_1,...,Kop_nZ], Kop_i an 'intop' functional
% declared by piesos_intopvar on the output function space of row i of Zx.
%
% Ground truths are independent of the function under test:
%  - decision variables: the number declared for row i is counted by hand
%    from the degrees: a functional of degree 1 in each of the k variables
%    on which row i depends, with cumulative degree d in each, has
%    (d+1)^k monomials, times the number of entries. Counts are read from
%    prog.var.idx, in the order of rows;
%  - values: for random values of the decision variables, p is evaluated
%    on random states x in two independent ways. Left: the operator
%    p.C.ops{1,j} acting on column j of Zx.C (a tensor monomial
%    x_k1 o ... o x_kd) is applied with apply_op to the product of the
%    factors, each renamed to its own copy of the variables, as 'intop'
%    expects. Right: the functionals Kop_i that piesos_polyvar should
%    declare (replayed by direct calls of piesos_intopvar with
%    hand-written arguments, on a copy of the same program so that the
%    decision variables coincide) are applied to the rows z_i of the basis
%    evaluated by apply_poly, as pointwise products of the states. This
%    tests the composition Cop*Zx, the restriction of the tensor monomials
%    to the diagonal, and the order of rows and columns;
%  - intmat and degrees: the formats of piesos_polyvar are compared by the
%    counts above, and with direct calls of piesos_intopvar, for rows of
%    output degree 1 and, for a basis built by hand, of output degree 2,
%    where the options 'full', 'nomult', 'lower' and 'upper' differ;
%  - SDP: p is used in actual programs (sossolve, SeDuMi), where
%    p(x_j)=F(x_j) for states x_j, and F a polynomial functional of known
%    kernels, is feasible exactly if p can represent F, with F recovered.
%
% Requires only that the functional of row i is declared by one call of
% piesos_intopvar.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_polyvar
%
% Copyright (C) 2026 PIETOOLS Team
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - The ground-truth value check depends on apply_op matching dummy
%   variables by name: monomial_value renames each factor's variables to
%   s_v_c (c = copy index) when the monomial has degree >=2 in s_v;
% - The Zd/kernels of intop are read through kernels(op) (dvar_names) and
%   Zd (rand_dvals); intop.params are not kernels;
% - case_sdp_recovery calls sossolve (SeDuMi); the equality-only program
%   gets a trivial PSD matrix because SeDuMi fails ("Size b mismatch")
%   without one.

tests ={@case_structure, @case_values, @case_degs_formats, @case_full_degs, ...
         @case_intmat_formats, @case_output_degree2, @case_mixed_basis, @case_three_variables_shape, @case_discarded_variable_terms, ...
         @case_three_variables, @case_vector_state, @case_constant_basis, ...
         @case_sequential, @case_piesos_eq, @case_sdp_recovery, @case_errors};
nfail = 0;
for i=1:numel(tests)
    try
        tests{i}();
        fprintf('PASS  %s\n',func2str(tests{i}));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_polyvar: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_polyvar: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function m = mx(maxdeg)
% Maximal degrees in the format of piesos_intopvar (braces: a cell must not
% expand to a struct array)
m = struct('max',{maxdeg});
end

function f = fl(c)
% Full degrees (monomials) in the format of piesos_intopvar
f = struct('full',{c});
end

function A = cum_degs(c,d)
% Degrees (rows) of all monomials of cumulative degree at most c in d variables
G = cell(1,d);
[G{:}] = ndgrid(0:c);
A = cell2mat(cellfun(@(g) g(:),G,'UniformOutput',false));
A = sortrows(A(sum(A,2)<=c,:));
end

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s".',err.message);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function [vartab,Zx,dep,vars] = basis_2d(degs)
% States x1(s1,s2) on [0,1]x[0,2], x2(s1) on [0,1], and the basis of monomials
% of cumulative degrees degs. dep(k,v): state k depends on s_v
pvar s1 s2
x1 = polyopvar('x1',[s1;s2],[0,1;0,2]);
x2 = polyopvar('x2',s1,[0,1]);
vartab = [x1;x2];
Zx = dmonomials(vartab,degs);
dep = [1 1; 1 0];
vars = {'s1','s2'};
end

function xvals = rand_states(dep,vars)
% Random polynomial states, depending on the variables dep(k,:)
xvals = cell(size(dep,1),1);
for k=1:size(dep,1)
    xvals{k} = rand_poly([1,1],vars(dep(k,:)>0),2,4);
end
end

function [K,prog] = replay(prog,Zx,intmat,degs)
% The functionals that piesos_polyvar should declare, one per row of Zx, by
% direct calls of piesos_intopvar with intmat{i}, degs{i} for the variables
% on which row i depends
n =size(Zx.C.ops,1);
K = cell(1,n);
for i=1:n
    is_dep = Zx.depmat.out(i,:)>0;
    [prog,K{i}] = piesos_intopvar(prog,[1,Zx.dims.out(i)],Zx.vars.out(is_dep),Zx.dom.out(is_dep,:),Zx.depmat.out(i,is_dep),intmat{i},degs{i});
end
end

function d0 = rand_dvals(Klist)
% Random values for the decision variables of the functionals in Klist (map by name)
d0 = containers.Map('KeyType','char','ValueType','double');
for j=1:numel(Klist)
    for l=1:numel(Klist{j}.Zd)
        d0(Klist{j}.Zd{l}) = randn;
    end
end
end

function v = monomial_value(Zrow,xvals,dep,vars)
% Value of the distributed monomial x_k1 o ... o x_kd (the nonzeros of Zrow)
% as the argument of a functional acting on it: the tensor product of the
% factors, where in factor f each variable s_v it depends on is renamed to the
% c-th copy of s_v, c = number of factors up to f depending on s_v. Copies are
% named s_v_c if the monomial has degree >=2 in s_v, and s_v otherwise
ks = full(Zrow);        ks = ks(ks>0);
tot = sum(dep(ks,:),1);
cnt = zeros(1,numel(vars));
v = polynomial(1);
for f=1:numel(ks)
    xk = xvals{ks(f)};
    for vv=find(dep(ks(f),:)>0)
        cnt(vv) = cnt(vv)+1;
        if tot(vv)>=2
            newname = [vars{vv},'_',num2str(cnt(vv))];
        else
            newname = vars{vv};
        end
        xk.varname(strcmp(xk.varname,vars{vv})) = {newname};
    end
    v = kron(v,xk);
end
end

function val = value_of_p(p,d0,Zx,xvals,dep,vars)
% p at the states xvals: sum over columns j of the operators in p, applied to
% the distributed monomial of column j
val = 0;
for j=1:size(p.C.ops,2)
    op = p.C.ops{1,j};
    if isempty(op)
        continue
    end
    opn = eval_intop_dpvar(op,d0);
    v = monomial_value(Zx.Z(j,:),xvals,dep,vars);
    val = val + double(apply_op(opn,v));
end
end

function val = value_of_K(K,d0,Zx,xvals)
% Sum over rows i of K{i} applied to row i of the basis at the states (pointwise products)
z = apply_poly(Zx,xvals);
blk = [0;cumsum(Zx.dims.out(:))];
val = 0;
for i=1:numel(K)
    Kn = eval_intop_dpvar(K{i},d0);
    val = val + double(apply_op(Kn,z(blk(i)+1:blk(i+1))));
end
end

function nm = dvar_names(p)
% Names of the decision variables in the kernels of p
nm = cell(0,1);
for j=1:size(p.C.ops,2)
    op = p.C.ops{1,j};
    if isempty(op)
        continue
    end
    opk = kernels(op);
    for k=1:numel(opk)
        if isa(opk{k},'dpvar')
            nm = [nm; opk{k}.dvarname(:)];
        end
    end
end
end

function cnt = block_sizes(prog)
% Number of decision variables in each of the variables declared in prog
cnt = diff(cell2mat(prog.var.idx));
cnt = cnt(:);
end

function check_p_structure(p,Zx)
% p is a scalar polyopvar on the monomials of Zx, with an intop of the degrees
% of that monomial in each column on which a row of Zx acts
assert(isa(p,'polyopvar') && isequal(size(p),[1,1]),'p must be a scalar polyopvar.');
assert(isequal(full(p.Z),full(Zx.Z)),'p must be defined by the monomials of Zx.');
assert(isequal(size(p.C.ops),[1,size(Zx.C.ops,2)]),'Wrong size of the coefficients of p.');
for j=1:size(Zx.C.ops,2)
    is_used = any(~cellfun(@isempty,Zx.C.ops(:,j)));
    op = p.C.ops{1,j};
    if ~is_used
        assert(isempty(op),'Column %d of Zx is not used, so p must not act on it.',j);
        continue
    end
    assert(isa(op,'intop'),'Column %d of p must be an ''intop''.',j);
    deps_j = Zx.C.depmat.in(j,:);
    assert(isequal(p.C.depmat.in(j,:),deps_j),'Column %d of p must act on a monomial of degrees %s.',j,mat2str(deps_j));
    assert(isequal(op.depmat.in,deps_j(deps_j>0)),'Operator %d of p has the wrong degrees.',j);
    assert(isequal(size(op),[1,Zx.C.dims.in(j)]),'Operator %d of p has the wrong size.',j);
end
initialize(p);          % must not error
end


%% ---- cases ------------------------------------------------------------
function case_structure()
[vartab,Zx] = basis_2d([0;1;2]);
prog0 = piesos_program(vartab);
[prog,p] = piesos_polyvar(prog0,Zx,struct('degs',2));
check_p_structure(p,Zx);
% one functional per row, with (2+1)^k monomials for the k variables of the row
R = Zx.depmat.out>0;
assert(isequal(R,logical([0 0;1 1;1 0;1 1;1 1;1 0])),'Unexpected basis; test assumes the order of dmonomials.');
assert(isequal(block_sizes(prog),3.^sum(R,2)),'Wrong number of decision variables: %s',mat2str(block_sizes(prog)'));
assert(prog.var.num==6,'One variable block per row of Zx.');
% default options: constant kernels
[prog1,p1] = piesos_polyvar(prog0,Zx);
check_p_structure(p1,Zx);
assert(isequal(block_sizes(prog1),ones(6,1)),'Default degrees must give constant kernels.');
end

function case_values()
rng(1);
[vartab,Zx,dep,vars] = basis_2d([0;1;2]);
prog0 = piesos_program(vartab);
[~,p] = piesos_polyvar(prog0,Zx,struct('degs',2));
n = size(Zx.C.ops,1);
K = replay(prog0,Zx,repmat({'full'},n,1),repmat({2},n,1));
d0 = rand_dvals(K);
for trial=1:3
    xvals = rand_states(dep,vars);
    lhs = value_of_p(p,d0,Zx,xvals,dep,vars);
    rhs = value_of_K(K,d0,Zx,xvals);
    assert(abs(lhs-rhs)<=1e-9*max(1,abs(rhs)),'p(x) differs from sum of the functionals on the monomials: %g vs %g.',lhs,rhs);
end
end

function case_degs_formats()
rng(2);
[vartab,Zx,dep,vars] = basis_2d([0;1;2]);
prog0 = piesos_program(vartab);
R = Zx.depmat.out>0;            % rows depend on [0 0;1 1;1 0;1 1;1 1;1 0]
nZ = 6;
% each: degrees as specified, and the resulting cumulative degree in each
% variable for each row (nZ x 2)
cases = {
    {[],                          zeros(nZ,2)}                              % default
    {3,                           3*ones(nZ,2)}                             % scalar
    {[1 4],                       repmat([1 4],nZ,1)}                       % per variable, for all monomials
    {(0:5)',                      repmat((0:5)',1,2)}                       % a degree for each monomial
    {mx(2),                       2*ones(nZ,2)}                             % struct, for all monomials
    {mx({1,4}),                   repmat([1 4],nZ,1)}                       % struct with a cell, for all monomials
    {{0;1;{2,7};{1,3};[2 1];2},   [0 0;1 1;2 0;1 3;2 1;2 2]}                % cell, mixed formats, N-wide cell for a row in s1 only
    {{0;mx(1);mx({2,7});mx({1,3});mx({2,1});mx(2)}, [0 0;1 1;2 0;1 3;2 1;2 2]}   % cell of structs
    };
for c=1:size(cases,1)
    if isempty(cases{c,1}{1})
        opts = struct();
    else
        opts = struct('degs',{cases{c,1}{1}});
    end
    [prog,p] = piesos_polyvar(prog0,Zx,opts);
    check_p_structure(p,Zx);
    expected = prod((cases{c,1}{2}+1).^R,2);
    assert(isequal(block_sizes(prog),expected),'Case %d: wrong number of decision variables: %s, expected %s.',c,mat2str(block_sizes(prog)'),mat2str(expected'));
end
% values, for per-variable degrees
degs_rows = {{0;1;{2,7};{1,3};{2,1};2}, {0;mx(1);mx({2,7});mx({1,3});mx({2,1});mx(2)}, {0;1;2;[1 3];[2 1];2}};
degs_for_K = {0;1;2;{1,3};{2,1};2};     % row 3 depends on s1 only
K = replay(prog0,Zx,repmat({'full'},nZ,1),degs_for_K);
d0 = rand_dvals(K);
xvals = rand_states(dep,vars);
for c=1:numel(degs_rows)
    [~,p] = piesos_polyvar(prog0,Zx,struct('degs',{degs_rows{c}}));
    assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for cell degrees %d.',c);
end
end

function case_full_degs()
% Monomials of the kernels specified with the field 'full' for rows depending
% on [0 0;1 1;1 0;1 1;1 1;1 0], intmat 'full', a single integral per variable
% (deps=1): a cell of the terms (here one) with a joint array for all dummy
% variables, a 1 x N cell with the monomials of each variable (no columns for
% a variable the row does not depend on), and one row as a single monomial
rng(12);
[vartab,Zx,dep,vars] = basis_2d([0;1;2]);
prog0 = piesos_program(vartab);
degs = {mx(0); fl({[1 0;0 2]}); fl({[0;3]}); fl({[0;1;2],[0;4]}); fl({[0 0;1 1;2 0]}); fl({[0;1],[]})};
[prog,p] = piesos_polyvar(prog0,Zx,struct('degs',{degs}));
check_p_structure(p,Zx);
assert(isequal(block_sizes(prog),[1;2;2;6;3;2]),'Wrong number of decision variables: %s',mat2str(block_sizes(prog)'));
degs_K = {mx(0); fl({[1 0;0 2]}); fl({[0;3]}); fl({[0;1;2],[0;4]}); fl({[0 0;1 1;2 0]}); fl({[0;1]})};
K = replay(prog0,Zx,repmat({'full'},6,1),degs_K);
d0 = rand_dvals(K);
xvals = rand_states(dep,vars);
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for full degrees.');
% a single row of the array is a single monomial: s1 s2^2 in the rows depending on both variables
[prog,~] = piesos_polyvar(prog0,Zx,struct('degs',{{mx(0); fl({[1 2]}); mx(0); fl({[1 2]}); fl({[1 2]}); mx(0)}}));
assert(isequal(block_sizes(prog),[1;1;1;1;1;1]),'A single row must be a single monomial: %s',mat2str(block_sizes(prog)'));
end

function case_intmat_formats()
rng(3);
[vartab,Zx,dep,vars] = basis_2d([0;1;2]);
prog0 = piesos_program(vartab);
nZ = 6;
ref_counts = 3.^sum(Zx.depmat.out>0,2);
% a char, also as string, for all monomials, 'no_mult' as 'nomult'
for opt = {'full','lower','upper','nomult','no_mult'}
    [prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',opt{1},'degs',2));
    check_p_structure(p,Zx);
    assert(isequal(block_sizes(prog),ref_counts),'Wrong number of decision variables for intmat %s.',opt{1});
end
[prog,~] = piesos_polyvar(prog0,Zx,struct('intmat',"lower",'degs',2));
assert(isequal(block_sizes(prog),ref_counts),'A string must work as a char.');
% a cell, per monomial: a char, or a 1 x N cell whose elements for variables the
% monomial does not depend on are ignored
I2 = {[0,2],[0,2]};
intmat = {{[0 5],[1 2 3]};      % row 1: constant, dependent on no variables
          I2;                   % row 2 (s1,s2)
          {[0,2],zeros(1,0)};   % row 3 (s1): nothing for s2
          'lower';              % row 4
          I2;                   % row 5
          {[0,2],[7 7]}};       % row 6 (s1): arbitrary element for s2
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',2));
check_p_structure(p,Zx);
assert(isequal(block_sizes(prog),ref_counts),'Wrong number of decision variables for a cell of intmats.');
% row 2 (monomial x1) must be the integral over the full domain of both variables
j = find(Zx.Z(:,1)==1 & Zx.Z(:,2)==0);
assert(isequal(p.C.ops{1,j}.intmat,{[0,2],[0,2]}),'Wrong intmat in the operator of row 2.');
K =replay(prog0,Zx,{'full';'full';{[0,2]};'lower';{[0,2],[0,2]};{[0,2]}},repmat({2},nZ,1));
d0 = rand_dvals(K);
xvals = rand_states(dep,vars);
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for cell intmat.');
end

function Zx = basis_degree2_output(dom)
% A basis with a single monomial x o x (identity I2 = I1 o I1), with output a
% function of two copies of s1, so that the functional acting on it has degree
% 2 and the intmat options differ
pvar s1
x = polyopvar('x',s1,dom);
I1 = x.C.ops{1,1};
I2 = otimes(I1,I1,[2,2]);
Zx = polyopvar();
Zx.C = tensopmat({I2});
Zx.Z = sparse([1 1]);
Zx.varname = {'x'};
Zx = initialize(Zx);
end

function Zx = basis_tensor_identity(vars,dom)
% A basis with a single monomial x o x, with output of degree 2 in each of the
% variables vars (tensor product identity of x(vars) on x o x)
x = polyopvar('x',polynomial(vars(:)),dom);
I1 = x.C.ops{1,1};
I2 = otimes(I1,I1,[2,2]);
Zx = polyopvar();
Zx.C = tensopmat({I2});
Zx.Z = sparse([1 1]);
Zx.varname = {'x'};
Zx = initialize(Zx);
end

function case_output_degree2()
rng(4);
Zx = basis_degree2_output([0,1]);
assert(isequal(Zx.depmat.out,2),'Test assumes the output of degree 2.');
vartab = polyopvar('x',polynomial({'s1'}),[0,1]);
prog0 = piesos_program(vartab);
% monomials per term, by hand: for degs=1, a term with two free dummy variables
% has 3 monomials, the multiplier (one free variable) 2
expected = struct('full',[2;3;3],'nomult',[3;3],'no_mult',[3;3],'lower',3,'upper',3);
opts = fieldnames(expected);
for o=1:numel(opts)
    [prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',opts{o},'degs',1));
    assert(isequal(block_sizes(prog),expected.(opts{o})),'Option %s: wrong number of decision variables: %s',opts{o},mat2str(block_sizes(prog)'));
    optname = strrep(opts{o},'no_mult','nomult');
    [~,K] = piesos_intopvar(prog0,[1,1],{'s1'},[0,1],2,optname,1);
    assert(isequal(p.C.ops{1,1}.intmat,K.intmat),'Option %s: wrong intmat.',opts{o});
    xv = {rand_poly([1,1],{'s1'},2,3)};
    d0 = rand_dvals({K});
    assert(abs(value_of_p(p,d0,Zx,xv,1,{'s1'})-value_of_K({K},d0,Zx,xv))<1e-9,'Option %s: values differ.',opts{o});
end
% default: 'full'
[prog,~] = piesos_polyvar(prog0,Zx,struct('degs',1));
assert(isequal(block_sizes(prog),[2;3;3]),'The default intmat must be ''full'': %s',mat2str(block_sizes(prog)'));
% explicit intmat, with degrees per term: [0,2,0,3] (s1<=s2, 3 monomials of
% cumulative degree 1) and [0,3,1,1] (s2=s1, 3 monomials of degree <=2 in the
% free variable). Options are nZ x 1 cells (here nZ=1): the element is a 1 x N
% cell of intmats, the monomials of the terms a cell of the terms (a vector,
% for N=1), in a struct with field 'full'
intmat_opt = {{[0,2,0,3; 0,3,1,1]}};
degs_term = fl({[0 0;1 0;0 1];[0 0;1 0;2 0]});
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat_opt},'degs',degs_term));
assert(isequal(block_sizes(prog),[3;3]),'Wrong number of decision variables for explicit intmat and degrees.');
K = replay(prog0,Zx,{{[0,2,0,3; 0,3,1,1]}},{degs_term});
d0 = rand_dvals(K);
xv = {rand_poly([1,1],{'s1'},2,3)};
assert(abs(value_of_p(p,d0,Zx,xv,1,{'s1'})-value_of_K(K,d0,Zx,xv))<1e-9,'Values differ for explicit intmat.');
% the same, as the element of an nZ x 1 cell
[prog,~] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat_opt},'degs',{{degs_term}}));
assert(isequal(block_sizes(prog),[3;3]),'Wrong number of decision variables for a cell with the full degrees.');
% a single row of an array is a single monomial: s_dum1^2 for the second term
% (s2=s1), and 2 monomials for the first (upper: coupled dummy variables)
degs_mons = fl({[0 0;1 1];[2 0]});
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat_opt},'degs',degs_mons));
assert(isequal(block_sizes(prog),[2;1]),'Monomials: wrong number of decision variables: %s',mat2str(block_sizes(prog)'));
K = replay(prog0,Zx,{{[0,2,0,3; 0,3,1,1]}},{degs_mons});
d0 = rand_dvals(K);
xv = {rand_poly([1,1],{'s1'},2,3)};
assert(abs(value_of_p(p,d0,Zx,xv,1,{'s1'})-value_of_K(K,d0,Zx,xv))<1e-9,'Values differ for explicit monomials.');
% maximal degrees for all terms: a scalar, or a 1 x N cell, also in a struct
% (cumulative degree 1: 3 monomials in the first term, 2 in the second)
for degs = {1, {1}, mx(1), mx({1})}
    [prog,~] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat_opt},'degs',degs));
    assert(isequal(block_sizes(prog),[3;2]),'Wrong number of decision variables for explicit intmat and one degree.');
end
% the cell of the terms must have an element for each term
assert_err(@() piesos_polyvar(prog0,Zx,struct('intmat',{intmat_opt},'degs',fl({[0 0];[0 0];[0 0]}))),'ntrms(1) x ... x ntrms(N) cell');
end

function case_mixed_basis()
% Two rows on different variables: x o x (degree 2 in s1 only) and y (degree 1
% in s2 only): N=2 but each row depends on one variable, so the specifications
% for the other variable must be dropped, also inside the degrees of each term
rng(9);
pvar s1 s2
Zx2 = basis_degree2_output([0,1]);
y = polyopvar('y',s2,[0,2]);
Zx = [Zx2;y];
assert(isequal(Zx.vars.out,{'s1','s2'}) && isequal(Zx.depmat.out,[2 0;0 1]),'Unexpected mixed basis.');
dep = [1 0;0 1];    vars = {'s1','s2'};
prog0 = piesos_program(Zx);
% 'lower' for both rows (3 monomials each: cumulative degree 1, resp. 2); N-wide
% degrees, with arbitrary values for the variable the row does not depend on
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat','lower','degs',{{{1,7};{9,2}}}));
assert(isequal(block_sizes(prog),[3;3]),'Wrong number of decision variables: %s',mat2str(block_sizes(prog)'));
K = replay(prog0,Zx,{'lower';'lower'},{{1};{2}});
d0 = rand_dvals(K);
xvals = {rand_poly([1,1],{'s1'},2,3);rand_poly([1,1],{'s2'},2,3)};
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for the mixed basis.');
% explicit intmat for row 1, N-wide with an empty element for s2, and its
% monomials as a cell of the terms (2 x 1); term 1 of cumulative degree 1,
% term 2 (s2=s1) up to 2
T = {[0 0;1 0;0 1];[0 0;1 0;2 0]};
intmat = {{[0,2,0,3; 0,3,1,1],zeros(1,0)};'full'};
degs = {fl(T);{9,2}};
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',{degs}));
assert(isequal(block_sizes(prog),[3;3;3]),'Wrong number of decision variables: %s',mat2str(block_sizes(prog)'));
K = replay(prog0,Zx,{{[0,2,0,3; 0,3,1,1]};'full'},{fl(T);{2}});
d0 = rand_dvals(K);
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for the mixed basis, explicit intmat.');
% the same monomials for all terms of row 1: a 1 x N cell, empty for s2; the tie
% s2=s1 discards the monomials in s_dum2
degs = {fl({T{1},[]});{9,2}};
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',{degs}));
assert(isequal(block_sizes(prog),[3;2;3]),'Wrong number of decision variables for a 1 x N cell of monomials: %s',mat2str(block_sizes(prog)'));
K = replay(prog0,Zx,{{[0,2,0,3; 0,3,1,1]};'full'},{fl({T{1}});{2}});
d0 = rand_dvals(K);
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for a 1 x N cell of monomials.');
% cumulative degrees [1 5] for (s1,s2) with an explicit intmat for row 1, used
% for each of its terms: 'upper' with 2 free dummies (3 monomials), s2=s1 tie
% with 1 free dummy (2 monomials); row 2: 6
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',[1 5]));
assert(isequal(block_sizes(prog),[3;2;6]),'Wrong number of decision variables for cumulative degrees and explicit intmat: %s',mat2str(block_sizes(prog)'));
K = replay(prog0,Zx,{{[0,2,0,3; 0,3,1,1]};'full'},{1;5});
d0 = rand_dvals(K);
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for cumulative degrees and explicit intmat.');
end

function case_discarded_variable_terms()
% A row depending on s2 only (deps=[0 2], N=2) with intmat 'nomult' has ntrms=[1 2]:
% its cell of terms is 1 x 2, like a 1 x N cell of per-variable monomials. Elements of
% the first have sum(deps)=2 columns, of the second deps(i), so none for s1
rng(13);
Zxx = basis_tensor_identity({'s2'},[0,1]);
Zx = [polyopvar('y',polynomial({'s1'}),[0,1]);Zxx];
assert(isequal(Zx.vars.out,{'s1','s2'}) && isequal(Zx.depmat.out,[1 0;0 2]),'Unexpected basis.');
assert(isequal(Zx.varname,{'x';'y'}),'Unexpected order of the states.');
dep = [0 1;1 0];    vars = {'s1','s2'};     % x(s2), y(s1)
prog0 = piesos_program(Zx);
T = {[0 0;1 0]; [0 0;0 1;1 1]};                % 2 and 3 monomials for the 2 terms of 'nomult'
% cell of terms for row 2, N-wide (1 x 2): each term gets its own monomials
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat','nomult','degs',{{0;fl(T.')}}));
assert(isequal(block_sizes(prog),[1;2;3]),'The cell of terms must not be read as a 1 x N cell: %s',mat2str(block_sizes(prog)'));
K = replay(prog0,Zx,{'nomult';'nomult'},{0;fl(T)});
d0 = rand_dvals(K);
xvals = {rand_poly([1,1],{'s2'},2,3);rand_poly([1,1],{'s1'},2,3)};
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for a cell of terms with a discarded variable.');
% a 1 x N cell of the monomials of each variable has no columns for s1; the same monomials are used for both terms
[prog,~] = piesos_polyvar(prog0,Zx,struct('intmat','nomult','degs',{{0;fl({[],[0 0;0 1;1 1]})}}));
assert(isequal(block_sizes(prog),[1;3;3]),'A 1 x N cell must give the same monomials for all terms: %s',mat2str(block_sizes(prog)'));
% a column for s1 is rejected, in the 1 x 2 cell of terms, and in the cells of a term
assert_err(@() piesos_polyvar(prog0,Zx,struct('intmat','nomult','degs',{{0;fl({5,[0 0;0 1;1 1]})}})),'a column for each of the 2 dummy variables');
assert_err(@() piesos_polyvar(prog0,Zx,struct('intmat','nomult','degs',{{0;fl({{5,[0 0;0 1]},{[],[0 0;0 1]}})}})),'deps(i) columns');
end

function case_three_variables_shape()
% N=3: row 1 is a tensor product identity in s1 and s3 (degree 2 in both), row 2
% depends on s2 only. The intmat and the monomials of each term of row 1 are
% specified for all 3 variables, with singleton dimension for s2 in the cell of terms
rng(11);
Zxx = basis_tensor_identity({'s1','s3'},[0,1;0,2]);
y = polyopvar('y',polynomial({'s2'}),[0,1]);
Zx = [Zxx;y];
assert(isequal(Zx.vars.out,{'s1','s2','s3'}) && isequal(Zx.depmat.out,[2 0 2;0 1 0]),'Unexpected basis.');
dep = [1 0 1;0 1 0];    vars = {'s1','s2','s3'};
prog0 = piesos_program(Zx);
I1 = [0,3,0,1; 0,2,0,3];                        % s1: 2 terms
I3 = [0,3,1,1; 0,3,0,1; 0,2,0,3];               % s3: 3 terms
degs_wide = cell(2,1,3);                        % terms of (s1,s2,s3)
degs_dep = cell(2,3);                           % terms of (s1,s3)
for i1=1:2
    for i3=1:3
        A1 = cum_degs((i1-1)+(i3-1),2);
        A3 = cum_degs((i1-1)+(i3-1),2);
        degs_dep{i1,i3} = {A1,A3};
        if mod(i3,2)==1
            degs_wide{i1,1,i3} = {A1,[],A3};    % monomials of each variable, none for s2
        else
            % all dummy variables together: those of s1 first, then those of s3
            degs_wide{i1,1,i3} = [repelem(A1,size(A3,1),1), repmat(A3,size(A1,1),1)];
        end
    end
end
intmat = {{I1,zeros(1,0),I3};'full'};
[prog,p] = piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',{{fl(degs_wide);3}}));
[K,progR] = replay(prog0,Zx,{{I1,I3};'full'},{fl(degs_dep);3});
assert(isequal(block_sizes(prog),block_sizes(progR)),'Wrong number of decision variables: %s vs %s',mat2str(block_sizes(prog)'),mat2str(block_sizes(progR)'));
assert(numel(block_sizes(prog))==7,'6 terms in the first row, 1 in the second.');
d0 = rand_dvals(K);
xvals = {rand_poly([1,1],{'s1','s3'},2,3);rand_poly([1,1],{'s2'},2,3)};
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<1e-9,'Values differ for 3 variables with singleton dimensions.');
% a cell of terms of the dependent variables only (2 x 3), or of the wrong shape, is rejected
assert_err(@() piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',{{fl(degs_dep);3}})),'ntrms(1) x ... x ntrms(N) cell');
assert_err(@() piesos_polyvar(prog0,Zx,struct('intmat',{intmat},'degs',{{fl(permute(degs_wide,[3 2 1]));3}})),'ntrms(1) x ... x ntrms(N) cell');
end

function case_three_variables()
rng(5);
pvar s1 s2 s3
x1 = polyopvar('x1',[s1;s2;s3],[0,1;0,2;0,1]);
x2 = polyopvar('x2',s2,[0,2]);
x3 = polyopvar('x3',[s1;s3],[0,1;0,1]);
vartab = [x1;x2;x3];
Zx = dmonomials(vartab,[1;2]);
dep = [1 1 1; 0 1 0; 1 0 1];
vars = {'s1','s2','s3'};
prog0 = piesos_program(vartab);
n = size(Zx.C.ops,1);
R = Zx.depmat.out>0;
D = repmat([1 2 0],n,1);
[prog,p] = piesos_polyvar(prog0,Zx,struct('degs',[1 2 0]));
check_p_structure(p,Zx);
assert(isequal(block_sizes(prog),prod((D+1).^R,2)),'Wrong number of decision variables.');
% values: degrees [1 2 0] for the variables each row depends on
degs_K = cell(n,1);
for i=1:n
    degs_K{i} = num2cell(D(i,R(i,:)));
end
K = replay(prog0,Zx,repmat({'full'},n,1),degs_K);
d0 = rand_dvals(K);
xvals = rand_states(dep,vars);
assert(abs(value_of_p(p,d0,Zx,xvals,dep,vars)-value_of_K(K,d0,Zx,xvals))<=1e-9,'Values differ for 3 variables.');
end

function case_vector_state()
% A basis with a vector-valued monomial: w (size 2) and w.*w (size 2, the
% tensor product w o w of size 4, restricted to its diagonal)
rng(6);
wdims = 0;
for attempt=1:50
    % size of the state is random in {1,2}
    [w,~,wdims,wvars] = rand_polyopvar(1,2,1);
    if wdims==2
        break
    end
end
assert(isequal(wdims,2),'Could not generate a state of size 2.');
Zx = [w;w.*w];
assert(isequal(Zx.dims.out(:)',[2 2]),'Test assumes two rows of size 2.');
prog0 = piesos_program(w);
[prog,p] = piesos_polyvar(prog0,Zx,struct('degs',2));
check_p_structure(p,Zx);
assert(isequal(block_sizes(prog),[6;6]),'Wrong number of decision variables.');
K = replay(prog0,Zx,{'full';'full'},{2;2});
d0 = rand_dvals(K);
xvals = {rand_poly([2,1],wvars,2,4)};
lhs = value_of_p(p,d0,Zx,xvals,1,wvars);
rhs = value_of_K(K,d0,Zx,xvals);
assert(abs(lhs-rhs)<=1e-9*max(1,abs(rhs)),'Values differ for a vector-valued state: %g vs %g.',lhs,rhs);
end

function case_constant_basis()
pvar s1
x = polyopvar('x',s1,[0,1]);
Zx = dmonomials(x,0);
prog0 = piesos_program(x);
[prog,p] = piesos_polyvar(prog0,Zx,struct('degs',3));
assert(isequal(size(p.C.ops),[1,1]) && isa(p.C.ops{1,1},'intop') && isempty(p.C.ops{1,1}.vars.in),'p must be a constant functional.');
assert(isequal(block_sizes(prog),1),'A constant has a single decision variable.');
end

function case_sequential()
rng(7);
[vartab,Zx] = basis_2d([1;2]);
prog0 = piesos_program(vartab);
[prog1,p1] = piesos_polyvar(prog0,Zx,struct('degs',1));
[prog2,p2] = piesos_polyvar(prog1,Zx,struct('degs',1));
assert(prog2.var.num==2*prog1.var.num,'The second variable must add its own decision variables.');
n1 = dvar_names(p1);    n2 = dvar_names(p2);
assert(~isempty(n1) && ~isempty(n2),'Polynomials must depend on decision variables.');
assert(isempty(intersect(n1,n2)),'Decision variables of different polynomials must be distinct.');
end

function case_piesos_eq()
% The polynomial can be constrained by piesos_eq: adds constraints
[vartab,Zx] = basis_2d([1;2]);
prog0 = piesos_program(vartab);
[prog,p] = piesos_polyvar(prog0,Zx,struct('degs',1));
prog2 = piesos_eq(prog,p);
assert(prog2.expr.num>prog.expr.num,'piesos_eq must add constraints.');
end

function case_sdp_recovery()
% p(x_j)=F(x_j) for 12 random states x_j, with F(x) = c0 + int k1 x + int k2 x^2,
% for the basis [1;x;x^2], is feasible only if the kernels of p are of degree
% 1, and then recovers F
rng(8);
a = 0.5;    b = 2;
c0 = 0.7;   k1 = @(s) 1+2*s;    k2 = @(s) 3-s;
F = @(xfun) c0 + integral(@(s) k1(s).*xfun(s),a,b,'AbsTol',1e-13,'RelTol',1e-13) ...
              + integral(@(s) k2(s).*xfun(s).^2,a,b,'AbsTol',1e-13,'RelTol',1e-13);
pvar s1
x = polyopvar('x',s1,[a,b]);
Zx = dmonomials(x,(0:2)');
dep = 1;    vars = {'s1'};
% states: random polynomials of degree <=2
J = 12;
xs = cell(J,1);     coefs = zeros(J,3);
for j=1:J
    coefs(j,:) = randn(1,3);
    xs{j} = polynomial(coefs(j,:)',(0:2)',{'s1'},[1,1]);
end
Fvals = zeros(J,1);
for j=1:J
    cj = coefs(j,:);
    Fvals(j) = F(@(s) cj(1)+cj(2)*s+cj(3)*s.^2);
end
[ok,prog,p] = solve_p_equals(Zx,x,struct('degs',1),xs,Fvals,dep,vars);
assert(ok,'The functional with kernels of degree 1 must be feasible.');
% the recovered functional equals F on a new state
ps = piesos_getsol(prog,p);
xnew = polynomial([0.3;-1.1;0.8],(0:2)',{'s1'},[1,1]);
val = value_of_p_numeric(ps,Zx,{xnew},dep,vars);
ref = F(@(s) 0.3-1.1*s+0.8*s.^2);
assert(abs(val-ref)<1e-5*max(1,abs(ref)),'Recovered functional differs from F on a new state: %g vs %g.',val,ref);
% constant kernels cannot represent k1, k2 non-constant
ok = solve_p_equals(Zx,x,struct('degs',0),xs,Fvals,dep,vars);
assert(~ok,'Constant kernels cannot represent F.');
end

function [ok,prog,p] = solve_p_equals(Zx,x,opts,xs,Fvals,dep,vars)
prog = piesos_program(x);
[prog,p] = piesos_polyvar(prog,Zx,opts);
% SeDuMi cannot solve an equality-only program: add a trivial PSD matrix
[prog,~] = sosposmatr(prog,1);
for j=1:numel(xs)
    val = 0;
    for c=1:size(p.C.ops,2)
        op = p.C.ops{1,c};
        if isempty(op)
            continue
        end
        val = val + apply_op(op,monomial_value(Zx.Z(c,:),xs(j),dep,vars));
    end
    prog = soseq(prog,val-Fvals(j));
end
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
ok = (info.pinf==0 && info.dinf==0);
end

function val = value_of_p_numeric(ps,Zx,xvals,dep,vars)
% p with numeric kernels evaluated at the states
val = 0;
for c=1:size(ps.C.ops,2)
    op = ps.C.ops{1,c};
    if isempty(op)
        continue
    end
    val = val + double(apply_op(op,monomial_value(Zx.Z(c,:),xvals,dep,vars)));
end
end

function case_errors()
[vartab,Zx] = basis_2d([0;1;2]);
prog = piesos_program(vartab);
assert_err(@() piesos_polyvar(prog,3),'polyopvar');
assert_err(@() piesos_polyvar(prog,Zx,3),'struct');
assert_err(@() piesos_polyvar(prog,Zx,struct('intmat',3)),'Option ''intmat''');
assert_err(@() piesos_polyvar(prog,Zx,struct('intmat',{{'full','full'}})),'Option ''intmat''');
assert_err(@() piesos_polyvar(prog,Zx,struct('intmat',{{1;1;1;1;1;1}})),'must be specified as ''char'' or 1 x N cell');
assert_err(@() piesos_polyvar(prog,Zx,struct('intmat','everything')),'Option for');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',[1 2 3 4])),'scalar, 1 x N array');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{1,2}})),'scalar, 1 x N array');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs','a')),'scalar, 1 x N array');
% elements for the wrong number of variables
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;[1 2 3];0;0;0;0}})),'Maximal degrees must be specified as a scalar, or as 1 x N cell');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;0;{1,2,3};0;0;0}})),'Maximal degrees must be specified as a scalar, or as 1 x N cell');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;{};0;0;0;0}})),'Maximal degrees must be specified as a scalar, or as 1 x N cell');
assert_err(@() piesos_polyvar(prog,Zx,struct('intmat',{{{[0,2],[0,2]};{[0,2]};'full';'full';'full';'full'}})),'monomial 2 must be specified as ''char'' or 1 x N cell');
assert_err(@() piesos_polyvar(prog,Zx,struct('intmat',{{cell(1,0);'full';'full';'full';'full';'full'}})),'monomial 1 must be specified as ''char'' or 1 x N cell');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{@sin;0;0;0;0;0}})),'Maximal degrees must be specified as a scalar, or as 1 x N cell');
% invalid structs are passed on to piesos_intopvar
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',struct('foo',1))),'single field ''max'' or ''full''');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',struct('max',1,'full',{{}}))),'single field ''max'' or ''full''');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',struct('full',{1,2}))),'scalar struct');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;fl({[0 0 0]});0;0;0;0}})),'a column for each of the 2 dummy variables');
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;fl(1);0;0;0;0}})),'Full degrees must be specified as ntrms(1) x ... x ntrms(N) cell');
% a cell is read as maximal degrees: a cell of the terms is not accepted there
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;{{1,1}};0;0;0;0}})),'Maximal degrees must be specified as a scalar, or as 1 x N cell');
% errors of piesos_intopvar are passed on
assert_err(@() piesos_polyvar(prog,Zx,struct('degs',{{0;-1;0;0;0;0}})),'nonnegative integers');
% a valid call after errors
[~,p] = piesos_polyvar(prog,Zx);
assert(isa(p,'polyopvar'),'Valid call must work after errors.');
end

