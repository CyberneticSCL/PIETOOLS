function Test_piesos_intopvar()
% TEST_PIESOS_INTOPVAR tests, using assertions, piesos_intopvar,
%   [prog,Kop] = piesos_intopvar(prog,matdim,vars,dom,deps,intmat,degs),
% which declares an 'intop' functional Kop whose kernels are matrices of
% polynomials in the dummy variables si_dumj, with decision variable
% coefficients.
%
% Ground truths are independent of the function under test:
%  - intmat: the rows generated for 'full', 'nomult', 'lower' and 'upper'
%    are compared with hand-written rows (deps<=2, and the chains for
%    deps=3), and the number of terms with the known counts (ordered Bell
%    numbers 1,3,13,75 for 'full', d! for 'nomult');
%  - regions: each term with kernel 1 is applied to a random polynomial v,
%    and compared with a nested Gauss-Legendre integral of v over the weak
%    ordering that the term should represent. The orderings are
%    enumerated by brute force here, so this tests the semantics of the
%    rows (including the deltas of ties) rather than their text;
%  - kernels: the monomials of each kernel are compared with a brute-force
%    enumeration of the exponents allowed by the degrees (cumulative,
%    per-variable maximal, or explicit; the same for all terms, or per
%    term), with monomials in dummy variables set equal to another one
%    discarded. Each decision variable must be the coefficient of exactly
%    one monomial in exactly one entry of one kernel (unit vectors, via
%    eval_dpvar_at);
%  - semantics: for random values of the decision variables, Kop is
%    applied (apply_op) to a random polynomial vector v, and compared with
%    the sum over all terms of the Gauss-Legendre integral of K*v over the
%    product of the weak orderings along the variables, including
%    unsorted variable names, non-unit domains, matrix-valued kernels, and
%    several variables;
%  - SDP: the declared Kop is used in actual programs (sossolve, SeDuMi),
%    where the equalities Kop*v_ij = F(v_ij) for a target functional F of
%    known kernel are feasible exactly if the kernels can represent F,
%    with the target kernel recovered from the solution.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_intopvar
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
% - intop.params are sparse coefficients, not kernels: kernels are read
%   via kernels(K) (helpers ker, kers, ndv) and decision variables via
%   K.Zd (rand_dvals);
% - apply_op matches dummy variables by name: names_of builds the names of
%   v (s_dumj, or the state variable if a direction has one dummy) in the
%   order of the dummy variables, so v must be defined in those names;
% - The oracle uses 8-node Gauss-Legendre quadrature, exact for the
%   polynomial integrands used here; rng(k) is fixed per case;
% - The SDP cases use sossolve (SeDuMi); equality-only programs get a
%   trivial PSD matrix (sosposmatr), as SeDuMi fails ("Size b mismatch")
%   without one.

tests ={@case_intmat_rows, @case_term_counts, @case_regions_1var, ...
         @case_monomials_formats, @case_monomials_per_term, ...
         @case_degs_term_formats, @case_degs_term_unsorted, @case_degs_term_shape, ...
         @case_ties_canonical, @case_generic_kernels, ...
         @case_semantics_1var, @case_semantics_2var_unsorted, ...
         @case_semantics_3var, @case_user_intmat, @case_no_variables, ...
         @case_struct_inputs, @case_sequential, @case_arithmetic, ...
         @case_sdp_recovery, ...
         @case_sdp_two_variables, @case_errors};
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
    error('Test_piesos_intopvar: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_intopvar: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function prog = new_prog(varnames)
% Fresh program in the given spatial variables
N = numel(varnames);
x1 = polyopvar('x1',polynomial(varnames(:)),[zeros(N,1),(1:N)']);
prog = piesos_program(x1);
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

function R = ranks_expected(d,opt)
% Rank vectors (rows) of the integrals that opt should include, for d dummy
% variables: R(t,j) is the rank of si_j in the increasing order of values,
% equal ranks meaning equal values. All d^d vectors are enumerated, and those
% using each of the ranks 1,...,max are kept, ordered as documented: fewest
% distinct values first, then decreasing lexicographic.
switch opt
    case 'lower'
        R = d:-1:1;
    case 'upper'
        R = 1:d;
    otherwise
        idx = (0:d^d-1)';
        R = zeros(d^d,d);
        for j=1:d
            R(:,j) = mod(floor(idx/d^(j-1)),d)+1;
        end
        keep = false(size(R,1),1);
        for t=1:size(R,1)
            keep(t) = isequal(unique(R(t,:)),1:max(R(t,:)));
        end
        R = R(keep,:);
        if strcmp(opt,'nomult')
            R = R(max(R,[],2)==d,:);
        end
        [~,o] = sortrows([max(R,[],2),-R]);
        R = R(o,:);
end
end

function cp = compile_poly(P,names)
% Prepare the m x n 'polynomial' (or double) P for evaluation at points x,
% where variable names{k} has the value x(k)
P = polynomial(P);
cp.m = size(P,1);       cp.n = size(P,2);
D = zeros(size(P.degmat,1),numel(names));
if ~isempty(P.varname)
    [tf,loc] = ismember(P.varname,names);
    assert(all(tf),'Polynomial depends on an unexpected variable.');
    D(:,loc) = full(P.degmat);
end
cp.D = D;
% 'polynomial' layout: monomial rows, columns (j-1)*m+i for entry (i,j)
cp.C = full(P.C).';
end

function M = eval_poly(cp,x)
mon = prod(x.^cp.D,2);
M = reshape(cp.C*mon,cp.m,cp.n);
end

function M = evalpoly(P,names,x)
M = eval_poly(compile_poly(P,names),x);
end

function v = rand_polyvec(names,n,nterms,maxdeg)
% n x 1 'polynomial' with random coefficients in the variables names, and
% nterms distinct random monomials of degree <=maxdeg in each variable
nv = numel(names);
E = unique(randi([0,maxdeg],nterms,nv),'rows');
v = polynomial(randn(size(E,1),n),E,names(:),[n,1]);
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

function Kn = eval_random(K)
% K with its 'dpvar' kernels evaluated at random decision variable values
Kn = eval_intop_dpvar(K,rand_dvals({K}));
end

function [dnames,vnames] = names_of(K)
% Names of the dummy variables defining the kernels, and of the variables of
% the polynomial v that apply_op expects, in the same order
dnames = cell(1,0);     vnames = cell(1,0);
for i=1:numel(K.vars.in)
    dnames = [dnames, K.varnames.in_dum{i}];
    if numel(K.varnames.in_dum{i})==1
        vnames = [vnames, K.vars.in(i)];
    else
        vnames = [vnames, K.varnames.in{i}];
    end
end
end

function val = chain_integral(f,Rs,lo,hi,n)
% Integral of f over the product of the chains lo(i)<=t_i(1)<=...<=t_i(Rs(i))
% <=hi(i), i=1,...,N, by nested Gauss-Legendre quadrature with n nodes; f takes
% a 1 x N cell, with the ith element the vector t_i.
[x,w] = gauss_legendre(n);
N = numel(Rs);
tcell = arrayfun(@(R) zeros(1,R),Rs,'UniformOutput',false);
if N==0
    val = f(tcell);
    return
end
val = nest(1,Rs(1),hi(1),tcell);
    function out = nest(i,k,up,tc)
        out = 0;
        tau = lo(i)+(up-lo(i))*(x+1)/2;
        wt = (up-lo(i))/2*w;
        for q=1:n
            tc{i}(k) = tau(q);
            if k>1
                out = out + wt(q)*nest(i,k-1,tau(q),tc);
            elseif i<N
                out = out + wt(q)*nest(i+1,Rs(i+1),hi(i+1),tc);
            else
                out = out + wt(q)*f(tc);
            end
        end
    end
end

function [x,w] = gauss_legendre(n)
% Nodes and weights on [-1,1] (Golub-Welsch)
k = 1:n-1;
beta = k./sqrt(4*k.^2-1);
[V,D] = eig(diag(beta,1)+diag(beta,-1));
[x,idx] = sort(diag(D));
w = 2*V(1,idx)'.^2;
end

function [tot,terms] = oracle_functional(Kn,v,dnames,vnames,Rexp,dom,matdim)
% Value of the functional Kn (numeric kernels) on v: the sum over all terms of
% the integral of K*v over the product of the orderings Rexp{i}(t_i,:), by
% nested quadrature, independently of apply_op
N = numel(Rexp);
deps = cellfun(@(a) size(a,2),Rexp);
ntrms = cellfun(@(a) size(a,1),Rexp);
Rmax = cell(1,N);
for i=1:N
    Rmax{i} = max(Rexp{i},[],2);
end
dpos = [0,cumsum(deps)];
vc = compile_poly(v,vnames);    % vnames match dnames by position
tot = zeros(matdim(1),1);
terms = cell([ntrms,1]);
sz = [ntrms,1];
idx = cell(1,numel(sz));
for k=1:numel(Kn.params)
    [idx{:}] = ind2sub(sz,k);
    r = cell(1,N);      Rs = zeros(1,N);
    for i=1:N
        r{i} = Rexp{i}(idx{i},:);
        Rs(i) = Rmax{i}(idx{i});
    end
    Kk = compile_poly(ker(Kn,k),dnames);
    f = @(tc) integrand(tc);
    val = chain_integral(f,Rs,dom(:,1)',dom(:,2)',8);
    terms{k} = val;
    tot = tot + val;
end
    function out = integrand(tc)
        s = zeros(1,dpos(end));
        for ii=1:N
            s(dpos(ii)+(1:deps(ii))) = tc{ii}(r{ii});
        end
        out = eval_poly(Kk,s)*eval_poly(vc,s);
    end
end

function check_semantics(prog,vars,dom,deps,opt,degs,matdim,Rexp_dirs,tol)
% Declare Kop by piesos_intopvar, evaluate its kernels at random decision
% variables, and compare apply_op with the quadrature oracle
[prog,K] = piesos_intopvar(prog,matdim,vars,dom,deps,opt,degs);
assert(isvalid(K),'Declared functional is not valid.');
% variables, domains, degrees must be those given, for the sorted variables
[vars_sorted,sidx] = sort(vars);
dom_sorted = dom(sidx,:);
assert(isequal(K.vars.in,vars_sorted),'Variables must be sorted.');
assert(isequal(K.dom.in,dom_sorted),'Each variable must keep its own domain.');
assert(isequal(K.depmat.in,reshape(deps(sidx),1,[])),'Each variable must keep its own degree.');
Rexp = Rexp_dirs(sidx);                 % orderings, for sorted variables
Kn = eval_random(K);
[dnames,vnames] = names_of(K);
v = rand_polyvec(vnames,matdim(2),5,2);
c = apply_op(Kn,v);
cnum = double(c);
cor = oracle_functional(Kn,v,dnames,vnames,Rexp,dom_sorted,matdim);
assert(norm(cnum(:)-cor(:))<=tol*max(1,norm(cor)),...
    'apply_op differs from the quadrature oracle: %g vs %g.',norm(cnum(:)),norm(cor(:)));
end


%% ---- cases ------------------------------------------------------------
function case_intmat_rows()
% Hand-written rows, in the order multiplier, lower, upper for deps=2
prog = new_prog({'s'});
expected = struct();
expected.full_1 = [0,2];
expected.full_2 = [0,3,1,1; 0,3,0,1; 0,2,0,3];
expected.nomult_2 = [0,3,0,1; 0,2,0,3];
expected.lower_2 = [0,3,0,1];
expected.upper_2 = [0,2,0,3];
expected.lower_3 = [0,4,0,1,0,2];       % s1>=s2>=s3: s1 over [a,b], s2 in [a,s1], s3 in [a,s2]
expected.upper_3 = [0,2,0,3,0,4];       % s1<=s2<=s3
expected.nomult_1 = [0,2];
expected.lower_1 = [0,2];
expected.upper_1 = [0,2];
fn = fieldnames(expected);
for k=1:numel(fn)
    parts = strsplit(fn{k},'_');
    opt = parts{1};     d = str2double(parts{2});
    [~,K] = piesos_intopvar(prog,[1,1],{'s'},[0,1],d,opt,0);
    assert(isequal(K.intmat{1},expected.(fn{k})),'Wrong intmat for %s, deps=%d.',opt,d);
end
% The ties in 'full' for deps=3, written by hand: variables set equal are
% taken equal to the first of them
[~,K3] = piesos_intopvar(prog,[1,1],{'s'},[0,1],3,'full',0);
assert(ismember([0,4,1,1,1,1],K3.intmat{1},'rows'),'Missing s1=s2=s3.');
assert(ismember([0,2,0,4,1,1],K3.intmat{1},'rows'),'Missing s1=s3<s2.');
assert(ismember([0,3,1,1,0,4],K3.intmat{1},'rows'),'Missing s1=s2<s3.');
assert(ismember([0,4,0,1,2,2],K3.intmat{1},'rows'),'Missing s3=s2<s1.');
end

function case_term_counts()
prog = new_prog({'s'});
fubini = [1,3,13,75];
for d=1:4
    [~,K] = piesos_intopvar(prog,[1,1],{'s'},[0,1],d,'full',0);
    assert(size(K.params,1)==fubini(d) && size(K.intmat{1},1)==fubini(d),'Wrong number of terms for ''full'', deps=%d.',d);
    assert(size(unique(K.intmat{1},'rows'),1)==fubini(d),'Duplicate rows for ''full'', deps=%d.',d);
    [~,K] = piesos_intopvar(prog,[1,1],{'s'},[0,1],d,'nomult',0);
    assert(size(K.params,1)==factorial(d),'Wrong number of terms for ''nomult'', deps=%d.',d);
end
% Terms along several variables form a cell with one dimension per variable
prog2 = new_prog({'s1','s2','s3'});
[~,K] = piesos_intopvar(prog2,[1,1],{'s1','s2','s3'},[0 1;0 1;0 1],[2,1,3],'full',0);
assert(isequal(size(K.params),[3,1,13]),'Wrong size of params.');
[~,K] = piesos_intopvar(prog2,[1,1],{'s1','s2','s3'},[0 1;0 1;0 1],[2,2,2],'nomult',0);
assert(isequal(size(K.params),[2,2,2]),'Wrong size of params.');
[~,K] = piesos_intopvar(prog2,[1,1],{'s1','s2','s3'},[0 1;0 1;0 1],[2,2,2],'lower',0);
assert(numel(K.params)==1,'''lower'' must give a single term.');
end

function case_regions_1var()
% Each term, with kernel 1, integrates v over the weak ordering it should
% represent: all options, deps=1,2,3, against quadrature
rng(2);
prog = new_prog({'s'});
a = 0.3;    b = 1.7;
opts = {'full','nomult','lower','upper'};
for d=1:3
    for o=1:numel(opts)
        [~,K] = piesos_intopvar(prog,[1,1],{'s'},[a,b],d,opts{o},0);
        R = ranks_expected(d,opts{o});
        assert(size(K.params,1)==size(R,1),'Number of terms differs for %s, deps=%d.',opts{o},d);
        [dnames,vnames] = names_of(K);
        v = rand_polyvec(vnames,1,6,2);
        for t=1:size(R,1)
            Kt_kernels = repmat({polynomial(0)},size(K.params));
            Kt_kernels{t} = polynomial(1);
            Kt = intop(Kt_kernels,K.intmat,K.vars,K.dom);
            cnum = double(apply_op(Kt,v));
            Rs = max(R(t,:));
            f = @(tc) evalpoly(v,vnames,tc{1}(R(t,:)));
            cor = chain_integral(f,Rs,a,b,8);
            assert(abs(cnum-cor)<=1e-10*max(1,abs(cor)),...
                'Term %d of %s, deps=%d: %g vs %g.',t,opts{o},d,cnum,cor);
        end
        % the terms must differ from one another
        if size(R,1)>1
            vals = zeros(size(R,1),1);
            for t=1:size(R,1)
                f = @(tc) evalpoly(v,vnames,tc{1}(R(t,:)));
                vals(t) = chain_integral(f,max(R(t,:)),a,b,8);
            end
            assert(min(pdist_min(vals))>1e-6,'Terms are not distinct.');
        end
    end
end
end

function dmin = pdist_min(vals)
% Smallest distance of each value to the others
n = numel(vals);
dmin = zeros(n,1);
for i=1:n
    dmin(i) = min(abs(vals(i)-vals([1:i-1,i+1:n])));
end
end

function monoms_ref = ref_monomials(deg_specs,is_sub,deps,mode)
% Brute-force the exponents (rows, over all dummy variables of all directions)
% of the monomials allowed by deg_specs{i}, where variables is_sub{i} are equal
% to another variable and so have exponent 0. Mode 'max': a spec is the maximal
% cumulative degree (scalar) or the maximal degree per variable (1 x d array);
% mode 'full': an array with a monomial in each row
N = numel(deps);
sets = cell(1,N);
for i=1:N
    d = deps(i);
    spec = deg_specs{i};
    is_max = strcmp(mode,'max');
    if is_max
        if isscalar(spec)
            mx = spec*ones(1,d);
        else
            mx = spec;
        end
    else
        mx = max(spec,[],1);
    end
    G = cell(1,d);
    vecs = arrayfun(@(m) 0:m,mx,'UniformOutput',false);
    [G{:}] = ndgrid(vecs{:});
    E = zeros(numel(G{1}),d);
    for j=1:d
        E(:,j) = G{j}(:);
    end
    if is_max
        if isscalar(spec)
            ok = sum(E,2)<=spec;
        else
            ok = all(E<=spec,2);
        end
    else
        ok = ismember(E,spec,'rows');
    end
    ok = ok & ~any(E(:,is_sub{i})>0,2);
    sets{i} = E(ok,:);
end
% Cartesian product
monoms_ref = zeros(1,0);
for i=1:N
    nA = size(monoms_ref,1);
    nB = size(sets{i},1);
    monoms_ref = [repelem(monoms_ref,nB,1),repmat(sets{i},nA,1)];
end
monoms_ref = sortrows(monoms_ref);
end

function check_kernel_monomials(K,degfun,deps,mode)
% Compare the monomials of each kernel of K with the brute-force enumeration;
% degfun(k,i) is the degree specification in direction i of term k, deps the
% degrees of the (sorted) variables, mode 'max' (default) or 'full' (see
% ref_monomials)
if nargin<4
    mode = 'max';
end
N = numel(K.vars.in);
[dnames,~] = names_of(K);
sz = [size(K.params,1:N),1];
idx = cell(1,numel(sz));
for k=1:numel(K.params)
    [idx{:}] = ind2sub(sz,k);
    specs = cell(1,N);      is_sub = cell(1,N);
    for i=1:N
        specs{i} = degfun(k,i);
        row = K.intmat{i}(idx{i},:);
        is_sub{i} = row(1:2:end)==row(2:2:end);
    end
    ref = ref_monomials(specs,is_sub,deps,mode);
    P = ker(K,k);
    if isempty(ref)
        assert(~isa(P,'dpvar') || isempty(P.dvarname),'Kernel must be zero.');
        continue
    end
    assert(isa(P,'dpvar'),'Kernel must be a dpvar.');
    [tf,loc] = ismember(P.varname,dnames);
    assert(all(tf),'Kernel %d depends on a variable that is not a dummy variable of Kop.',k);
    got = zeros(size(P.degmat,1),numel(dnames));
    got(:,loc) = full(P.degmat);
    got = sortrows(got);
    assert(isequal(got,ref),'Monomials of kernel %d differ from the reference.',k);
end
end

function case_monomials_formats()
rng(3);
prog = new_prog({'s1','s2'});
vars = {'s1','s2'};
dom = [0 1;0 1];
% (d1,d2) = (2,1), 'full' (generated intmat): 3 and 1 terms. A cell or array
% is read as maximal degrees: struct('max',degs)
cases = {{0,0},{2,1},{[1 2],3},{[2 2],[1]}};
for c=1:numel(cases)
    degs = cases{c};
    [~,K] = piesos_intopvar(prog,[2,3],vars,dom,[2,1],'full',degs);
    check_kernel_monomials(K,@(k,i) degs{i},[2,1],'max');
    [~,K2] = piesos_intopvar(prog,[2,3],vars,dom,[2,1],'full',dmax(degs));
    check_kernel_monomials(K2,@(k,i) degs{i},[2,1],'max');
end
% a scalar: maximal cumulative degree for all variables and terms
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,2],'nomult',2);
check_kernel_monomials(K,@(k,i) 2,[2,2],'max');
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,2],'nomult',dmax(2));
check_kernel_monomials(K,@(k,i) 2,[2,2],'max');
% Default degrees: constant kernels
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],'full');
check_kernel_monomials(K,@(k,i) 0,[2,1],'max');
% monomials as 1 x N cell, per variable (2 and 1 columns), for each term of a generated intmat
dm = {[0 0;1 1;2 0],[0;2;3]};
[~,K] = piesos_intopvar(prog,[2,3],vars,dom,[2,1],'full',dfull(dm));
check_kernel_monomials(K,@(k,i) dm{i},[2,1],'full');
% a single row is a single monomial s1_dum1 s1_dum2^2 s2_dum1^3, not maximal degrees
% (kernel 1 has a tie between s1_dum1 and s1_dum2: zero; kernel 2 is lower)
d1 = {[1 2],3};
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],'full',dfull(d1));
check_kernel_monomials(K,@(k,i) d1{i},[2,1],'full');
assert(isequal(monomials_of(K,2),[1 2 3]) && ndv(K,2)==1,'A single row must be a single monomial, not maximal degrees.');
% repeated rows are a single monomial
[~,K] = piesos_intopvar(prog,[1,1],{'s1'},[0 1],1,'full',dfull({[0;1;1;2]}));
assert(ndv(K,1)==3,'Repeated rows must give a single monomial: %d decision variables.',ndv(K,1));
% an explicit intmat
intmat = {[0,3,0,1; 0,2,0,3],[0,2]};
ds = {[1 1],2};
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,ds);
check_kernel_monomials(K,@(k,i) ds{i},[2,1],'max');
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,1);
check_kernel_monomials(K,@(k,i) 1,[2,1],'max');
% monomials of each term (2 terms along s1, 1 along s2), as 1 x N cell
dm2 = {[0 0;1 1],[0;2]};
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,dfull({dm2;dm2}));
check_kernel_monomials(K,@(k,i) dm2{i},[2,1],'full');
assert(ndv(K,1)==4,'Products of 2 and 2 monomials give 4 monomials.');
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,dfull({d1;d1}));
assert(isequal(monomials_of(K,1),[1 2 3]) && ndv(K,1)==1,'A single row must be a single monomial, not maximal degrees.');
% monomials cannot be given as maximal degrees
assert_err(@() piesos_intopvar(prog,[2,3],vars,dom,[2,1],'full',{[0 1;1 1;2 0;0 2],[0;2;3]}),'Use the field ''full''');
assert_err(@() piesos_intopvar(prog,[2,3],vars,dom,[2,1],intmat,dmax({[0 1;1 1;2 0;0 2],[0;2;3]})),'Use the field ''full''');
end

function case_monomials_per_term()
% Monomials specified per term, with unsorted variables. The call has variables
% (s2,s1), with 2 terms along s2 and 3 along s1: a 2 x 3 cell, each element a
% cell with arrays for (s2,s1). The functional has variables (s1,s2), so params
% is 3 x 2, and the monomials of term (i1,i2) along (s1,s2) must be those given
% for term (i2,i1) along (s2,s1)
rng(4);
prog = new_prog({'s1','s2'});
vars = {'s2','s1'};
dom = [0 1;0 2];                   % s2 on [0,1], s1 on [0,2]
intmat = {[0,3,0,1; 0,2,0,3],[0,3,1,1; 0,3,0,1; 0,2,0,3]};     % for s2 (2 terms), s1 (3 terms)
degs_sorted = cell(3,2);           % degs_sorted{i1,i2}{i}: monomials along (s1,s2)
for k=1:numel(degs_sorted)
    degs_sorted{k} = {unique(randi([0,2],2,2),'rows'),unique(randi([0,2],2,2),'rows')};
end
degs_user = cell(2,3);
for i1=1:3
    for i2=1:2
        degs_user{i2,i1} = {degs_sorted{i1,i2}{2},degs_sorted{i1,i2}{1}};
    end
end
[~,K] = piesos_intopvar(prog,[1,2],vars,dom,[2,2],intmat,dfull(degs_user));
assert(isequal(K.vars.in,{'s1','s2'}),'Variables must be sorted.');
assert(isequal(K.dom.in,[0 2;0 1]),'Domains must follow the sorted variables.');
assert(isequal(size(K.params),[3,2]),'Wrong size of params.');
assert(isequal(K.intmat{1},intmat{2}) && isequal(K.intmat{2},intmat{1}),'Intmats must follow the sorted variables.');
check_kernel_monomials(K,@(k,i) degs_sorted{k}{i},[2,2],'full');
end

function case_ties_canonical()
% A kernel cannot depend on a dummy variable set equal to another one: for
% s2=s1 the monomials in s_dum2 are discarded. By hand for deps=2, cumulative
% degree 2: 6 monomials in (s_dum1,s_dum2) per integral, 3 in s_dum1 for the multiplier
prog = new_prog({'s'});
[~,K] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',2);
assert(isequal(K.intmat{1}(1,:),[0,3,1,1]),'Expected the multiplier first.');
nd = zeros(1,3);
for k=1:3
    nd(k) = ndv(K,k);
end
assert(isequal(nd,[3,6,6]),'Wrong number of decision variables: %s',mat2str(nd));
assert(~ismember('s_dum2',ker(K,1).varname),'Multiplier kernel must not depend on s_dum2.');
% deps=3, cumulative degree 1: with f free variables, 1+f monomials
[~,K3] = piesos_intopvar(prog,[1,1],{'s'},[0,1],3,'full',1);
for k=1:size(K3.params,1)
    row = K3.intmat{1}(k,:);
    nfree = nnz(row(1:2:end)~=row(2:2:end));
    assert(ndv(K3,k)==1+nfree,'Term %d: expected %d monomials.',k,1+nfree);
end
% explicit monomials all of nonzero degree in s_dum2: zero multiplier kernel, no decision variables
nvar0 = prog.var.num;
mons_z = [0 1; 0 2];                            % monomials s_dum2, s_dum2^2
[progz,Kz] = piesos_intopvar(prog,[1,2],{'s'},[0,1],2,'full',dfull({mons_z;mons_z;mons_z}));  % multiplier, lower, upper
assert(isequal(size(Kz.params),[3,1]) && isequal(Kz.intmat{1}(1,:),[0,3,1,1]),'Unexpected terms.');
assert(isequal(size(ker(Kz,1)),[1,2]) && ~isa(ker(Kz,1),'dpvar'),'Multiplier kernel must be a zero polynomial.');
assert(~any(any(full(ker(Kz,1).C))),'Multiplier kernel must be zero.');
assert(isa(ker(Kz,2),'dpvar') && isa(ker(Kz,3),'dpvar'),'Integral kernels must contain decision variables.');
assert(progz.var.num==nvar0+2,'Only the two integral kernels may declare decision variables.');
assert(isvalid(Kz),'Functional with a zero kernel must be valid.');
end

function case_generic_kernels()
% Each decision variable is the coefficient of exactly one monomial in exactly
% one entry of one kernel: kernels evaluated at unit vectors give a single 1
prog = new_prog({'s1','s2'});
[prog,K] = piesos_intopvar(prog,[2,3],{'s1','s2'},[0 1;0 1],[2,1],'full',{[1 1],2});
nd_tot = 0;
names_all = {};
for k=1:numel(K.params)
    P = ker(K,k);
    nZ = size(P.degmat,1);
    assert(numel(P.dvarname)==2*3*nZ,'Kernel %d: expected m*n*nmons decision variables.',k);
    nd_tot = nd_tot+numel(P.dvarname);
    names_all = [names_all; P.dvarname(:)];
    Ccoef = zeros(2*3*nZ,2*3*nZ);
    for l=1:numel(P.dvarname)
        d0 = containers.Map(P.dvarname,num2cell(zeros(numel(P.dvarname),1)));
        d0(P.dvarname{l}) = 1;
        Pl = eval_dpvar_at(P,d0);
        Ccoef(:,l) = reshape(full(Pl.C),[],1);
    end
    assert(isequal(sum(Ccoef~=0,1),ones(1,2*3*nZ)),'Each decision variable must give one coefficient.');
    assert(isequal(sum(Ccoef~=0,2),ones(2*3*nZ,1)),'Each coefficient must come from one decision variable.');
    assert(all(Ccoef(Ccoef~=0)==1),'Coefficients must be 1.');
end
assert(numel(unique(names_all))==nd_tot,'Decision variables must be distinct across kernels.');
assert(prog.var.num==numel(K.params),'One decision variable block per kernel is expected.');
end

function case_semantics_1var()
rng(5);
opts = {'full','nomult','lower','upper'};
for d=1:3
    for o=1:numel(opts)
        prog = new_prog({'s'});
        R = {ranks_expected(d,opts{o})};
        check_semantics(prog,{'s'},[-0.5,1.5],d,opts{o},{1},[2,3],R,1e-9);
    end
end
% different degrees, zero-degree kernels
prog = new_prog({'s'});
check_semantics(prog,{'s'},[0,2],2,'full',{[0 2]},[1,1],{ranks_expected(2,'full')},1e-9);
end

function case_semantics_2var_unsorted()
rng(6);
% variables given as (s2,s1), different degrees and domains
prog = new_prog({'s1','s2'});
vars = {'s2','s1'};
dom = [0.2 1.1; -1 0.5];                            % s2 on [0.2,1.1], s1 on [-1,0.5]
Rexp_dirs = {ranks_expected(1,'full'),ranks_expected(2,'full')};      % deps (s2,s1) = (1,2)
check_semantics(prog,vars,dom,[1,2],'full',{1,1},[2,2],Rexp_dirs,1e-9);
Rexp_dirs = {ranks_expected(2,'nomult'),ranks_expected(2,'nomult')};   % (s2,s1)
check_semantics(prog,vars,dom,[2,2],'nomult',{1,[1 0]},[1,2],Rexp_dirs,1e-9);
end

function case_semantics_3var()
rng(7);
% three variables; the 'lower' option (one term along each) and 'full'
prog = new_prog({'s1','s2','s3'});
vars = {'s1','s2','s3'};
dom = [0 1; 0.5 1.5; -1 1];
Rexp_dirs = {ranks_expected(1,'lower'),ranks_expected(2,'lower'),ranks_expected(2,'lower')};
check_semantics(prog,vars,dom,[1,2,2],'lower',{1,1,0},[1,1],Rexp_dirs,1e-9);
Rexp_dirs = {ranks_expected(1,'full'),ranks_expected(2,'full'),ranks_expected(1,'full')};
check_semantics(prog,vars,dom,[1,2,1],'full',{0,0,1},[1,2],Rexp_dirs,1e-9);
end

function case_user_intmat()
% A user-specified intmat, with the rows written by hand together with the
% ordering they stand for ([0,2,0,3] is s1<=s2, rank (1,2); [0,3,1,1] is
% s2=s1, rank (1,1))
rng(8);
prog = new_prog({'s'});
intmat = {[0,2,0,3; 0,3,1,1]};
R = {[1 2;1 1]};
check_semantics(prog,{'s'},[0,1.3],2,intmat,{[1 2]},[2,1],R,1e-9);
% per-term degrees, same number of monomials
degs = {[0 0];[2 0;1 0]};
[~,K] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,intmat,dfull(degs));
check_kernel_monomials(K,@(k,i) degs{k},2,'full');
end

function m = dmax(v)
% Maximal degrees v in the format of piesos_intopvar (braces: v can be a cell)
m = struct('max',{v});
end

function m = dfull(v)
% Full degrees v in the format of piesos_intopvar (braces: v is a cell)
m = struct('full',{v});
end

function E = simplex_rows(D,d)
% Degrees (rows) of all monomials in D variables of cumulative degree at most d
G = cell(1,D);
[G{:}] = ndgrid(0:d);
E = zeros(numel(G{1}),D);
for j=1:D
    E(:,j) = G{j}(:);
end
E = E(sum(E,2)<=d,:);
end

function mons = monomials_of(K,k)
% Exponents (rows, sorted) of the monomials of kernel k of K, over all dummy
% variables of K in the order of K.varnames.in_dum
[dnames,~] = names_of(K);
P = ker(K,k);
[tf,loc] = ismember(P.varname,dnames);
assert(all(tf),'Kernel %d depends on a variable that is not a dummy variable of Kop.',k);
mons = zeros(size(P.degmat,1),numel(dnames));
mons(:,loc) = full(P.degmat);
mons = sortrows(mons);
end

function case_degs_term_formats()
% 'full': the monomials of each term, as ntrms(1) x ... x ntrms(N) cell, each
% element an array over all dummy variables, or a 1 x N cell of arrays per
% variable. Variables: s1 (deps 2, 3 terms: tie, lower, upper) and s2 (deps 1,
% 1 term); dummy order (s1_dum1,s1_dum2,s2_dum1). Monomials by hand
prog = new_prog({'s1','s2'});
vars = {'s1','s2'};     dom = [0 1;0 1];
intmat = {[0,3,1,1; 0,3,0,1; 0,2,0,3],[0,2]};
degs = cell(3,1);
degs{1} = [0 0 0; 0 0 1; 1 0 0; 1 0 1; 2 0 0; 2 0 1; 0 3 0];    % the last has s1_dum2, set equal to s1_dum1: discarded
degs{2} = [0 0 0; 1 1 1; 0 2 1];                                % coupling the variables
degs{3} = simplex_rows(3,1);                                    % cumulative degree at most 1 in all 3
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,dfull(degs));
assert(isequal(size(K.params),[3,1]),'Wrong size of params.');
assert(isequal(monomials_of(K,1),[0 0 0; 0 0 1; 1 0 0; 1 0 1; 2 0 0; 2 0 1]),'Term 1 (tie): wrong monomials.');
assert(isequal(monomials_of(K,2),[0 0 0; 0 2 1; 1 1 1]),'Term 2 (coupled monomials): wrong monomials.');
assert(isequal(monomials_of(K,3),sortrows(simplex_rows(3,1))),'Term 3: wrong monomials.');
% generated intmat, same 3 x 1 terms
[~,Kg] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],'full',dfull(degs));
for k=1:3
    assert(isequal(monomials_of(Kg,k),monomials_of(K,k)),'Cell of terms for a generated intmat: wrong monomials in term %d.',k);
end
% a single row is a single monomial; with the tie s1_dum2=s1_dum1 a monomial in
% s1_dum2 is discarded, leaving a zero kernel
degs1 = {[1 0 1]; [1 2 1]; [0 1 0]};
[~,K1] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,dfull(degs1));
assert(isequal(monomials_of(K1,1),[1 0 1]) && ndv(K1,1)==1,'Term 1: a single row must be a single monomial.');
assert(isequal(monomials_of(K1,2),[1 2 1]) && ndv(K1,2)==1,'Term 2: a single row must be a single monomial.');
assert(isequal(monomials_of(K1,3),[0 1 0]) && ndv(K1,3)==1,'Term 3: a single row must be a single monomial.');
[~,K1z] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,dfull({[1 2 1]; [1 2 1]; [1 2 1]}));
assert(~isa(ker(K1z,1),'dpvar'),'A monomial in s1_dum2, equal to s1_dum1, must give a zero kernel.');
% each element as 1 x N cell, with the monomials of s1 (2 columns) and s2 (1)
degs2 = {{[2 0;1 0],[0;1]}; {[0 0;1 1;0 2],[0;1]}; {[0 0],1}};
[~,K2] = piesos_intopvar(prog,[1,1],vars,dom,[2,1],intmat,dfull(degs2));
assert(isequal(monomials_of(K2,1),[1 0 0; 1 0 1; 2 0 0; 2 0 1]),'Term 1 as cell: wrong monomials.');
assert(isequal(monomials_of(K2,2),sortrows([0 0 0;0 0 1;1 1 0;1 1 1;0 2 0;0 2 1])),'Term 2 as cell: wrong monomials.');
assert(isequal(monomials_of(K2,3),[0 0 1]),'Term 3 as cell: a single monomial.');
% degs{i1,i2} is the monomials of term i1 along s1 and i2 along s2; intmat{2}: (lower,upper)
intmat2 = {[0,3,1,1; 0,3,0,1; 0,2,0,3],[0,3,0,1; 0,2,0,3]};
degs3 = cell(3,2);
for i1=1:3
    for i2=1:2
        degs3{i1,i2} = simplex_rows(4,i1+i2-2);     % cumulative degree i1+i2-2 in the 4 dummy variables
    end
end
[~,K3] = piesos_intopvar(prog,[1,1],vars,dom,[2,2],intmat2,dfull(degs3));
assert(isequal(size(K3.params),[3,2]),'Wrong size of params.');
% monomials of term (i1,i2): cumulative degree d=i1+i2-2 in the free dummies
% (the tie in s1, i1=1, leaves 1 free dummy for s1; s2 has 2)
nfree1 = [1,2,2];   nfree2 = [2,2];
for i1=1:3
    for i2=1:2
        d = i1+i2-2;
        nf = nfree1(i1)+nfree2(i2);
        expected = nchoosek(d+nf,nf);
        assert(size(monomials_of(K3,sub2ind([3,2],i1,i2)),1)==expected,'Term (%d,%d): expected %d monomials.',i1,i2,expected);
    end
end
% only term (1,2) (1 along s1, 2 along s2) gets degree 3
degs4 = repmat({zeros(1,4)},3,2);
degs4{1,2} = simplex_rows(4,3);
[~,K4] = piesos_intopvar(prog,[1,1],vars,dom,[2,2],intmat2,dfull(degs4));
nmons = cellfun(@(P) numel(P.dvarname),kers(K4));
expected = ones(3,2);
expected(1,2) = nchoosek(3+3,3);                  % 3 free dummies
assert(isequal(nmons,expected),'Only term (1,2) may have degree 3: %s',mat2str(nmons));
% for N=1, a vector (row or column) of length ntrms
prog1 = new_prog({'s'});
dd = {simplex_rows(2,0),simplex_rows(2,1),simplex_rows(2,2)};
for dv = {dd(:),dd(:)'}
    [~,Kv] = piesos_intopvar(prog1,[1,1],{'s'},[0,1],2,{[0,3,1,1; 0,3,0,1; 0,2,0,3]},dfull(dv{1}));
    Kvk = kers(Kv);
    nm = cellfun(@(P) numel(P.dvarname),Kvk(:))';
    % tie (1 free dummy, degree 0): 1; lower (2 free dummies, degree 1): 3; upper (degree 2): 6
    assert(isequal(nm,[1,3,6]),'Cumulative degrees 0,1,2: wrong number of monomials: %s',mat2str(nm));
end
end

function case_degs_term_unsorted()
% Variables given as (s2,s1), with monomials of a term as array over the dummy
% variables in the specified order: (s2_dum1,s1_dum1,s1_dum2)
prog = new_prog({'s1','s2'});
vars = {'s2','s1'};     dom = [0 2;0 1];
intmat = {[0,2],[0,3,0,1]};         % s2: 1 term; s1: lower
degs = {[0 0 0; 1 0 2; 0 1 1]};     % (s2_dum1,s1_dum1,s1_dum2)
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat,dfull(degs));
assert(isequal(K.vars.in,{'s1','s2'}) && isequal(K.dom.in,[0 1;0 2]),'Variables must be sorted.');
% in the order (s1_dum1,s1_dum2,s2_dum1): [0 0 0], [0 2 1], [1 1 0]
assert(isequal(monomials_of(K,1),[0 0 0; 0 2 1; 1 1 0]),'Columns of the degrees must follow the specified order of the variables.');
% a single row: a single monomial s2_dum1 s1_dum2^2, in the order (s1_dum1,s1_dum2,s2_dum1)
[~,K1] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat,dfull({[1 0 2]}));
assert(isequal(monomials_of(K1,1),[0 2 1]),'A single row must be a single monomial, with the columns reordered.');
% maximal degrees: for s2 cumulative 2, for s1 (1,0)
[~,K2] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat,{2,[1 0]});
assert(isequal(monomials_of(K2,1),sortrows([0 0 0;0 0 1;0 0 2;1 0 0;1 0 1;1 0 2])),'Maximal degrees as cell, unsorted variables: wrong monomials.');
% 1 x N cell of the monomials of each variable: s2 (1 column), s1 (2 columns)
[~,K3] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat,dfull({[0;1],[1 0;0 1]}));
assert(isequal(monomials_of(K3,1),[0 1 0;0 1 1;1 0 0;1 0 1]),'Monomials of each variable, unsorted variables: wrong monomials.');
end

function case_degs_term_shape()
% The shape of the cell of monomials per term must match ntrms; a 1 x N cell of
% per-variable monomials is told apart from a cell of terms by the number of
% columns of its elements
prog = new_prog({'s1','s2'});
vars = {'s1','s2'};     dom = [0 1;0 1];
intmat = {[0,3,0,1; 0,2,0,3],[0,3,1,1; 0,3,0,1; 0,2,0,3]};     % 2 x 3 terms
d23 = repmat({zeros(1,4)},2,3);
[~,K] = piesos_intopvar(prog,[1,1],vars,dom,[2,2],intmat,dfull(d23));
assert(isequal(size(K.params),[2,3]),'Wrong size of params.');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[2,2],intmat,dfull(d23')),'ntrms(1) x ... x ntrms(N) cell');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[2,2],intmat,dfull(repmat({zeros(1,4)},1,6))),'ntrms(1) x ... x ntrms(N) cell');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[2,2],intmat,dfull(repmat({zeros(1,4)},6,1))),'ntrms(1) x ... x ntrms(N) cell');
% N=2, deps=[1,2] and ntrms=[1,2]: a 1 x 2 cell is either of the two variables, or of the two terms
intmat12 = {[0,2],[0,3,0,1; 0,2,0,3]};
% monomials of s1 (1 column) and s2 (2 columns), for both terms
[~,Kv] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({[0;1;2],[0 0;1 1]}));
assert(isequal(cellfun(@(P) numel(P.dvarname),kers(Kv)),[6,6]),'Monomials of each variable: 3x2 monomials in each term.');
assert(isequal(monomials_of(Kv,1),[0 0 0;0 1 1;1 0 0;1 1 1;2 0 0;2 1 1]) && isequal(monomials_of(Kv,2),monomials_of(Kv,1)),'Monomials of each variable: wrong monomials.');
% the two terms, as arrays with 3 columns (single monomials)
[~,Kt] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({[2 1 2],[0 1 1]}));
assert(isequal(cellfun(@(P) numel(P.dvarname),kers(Kt)),[1,1]),'A single row must be a single monomial.');
assert(isequal(monomials_of(Kt,1),[2 1 2]) && isequal(monomials_of(Kt,2),[0 1 1]),'Wrong monomials of the terms.');
% the two terms, as 1 x N cells
[~,Kt2] = piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({{[0;1],[0 0;1 1]},{2,[1 2]}}));
assert(isequal(cellfun(@(P) numel(P.dvarname),kers(Kt2)),[4,1]),'Cells for each of the terms: wrong number of monomials.');
% N=1: a cell of the single term, or of the monomials of the variable, are the same
prog1 = new_prog({'s'});
[~,K1a] = piesos_intopvar(prog1,[1,1],{'s'},[0,1],2,{[0,3,0,1]},dfull({[0 0;1 0]}));
assert(ndv(K1a,1)==2,'One variable, one term.');
% errors: the width of the arrays matches neither
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({[2 1],[0 1 1]})),'column for each of the 3 dummy variables');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({[0;1],[0 1 1]})),'column for each of the 3 dummy variables');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({2,[0 1 1]})),'column for each of the 3 dummy variables');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({'a',[0 1 1]})),'array, or as 1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({{1},[0 1 1]})),'1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({[2 -1 2],[0 1 1]})),'nonnegative integers');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull({[0 1 1]})),'with an element for each term');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,dfull([0 1 1])),'Full degrees must be specified');
% the struct: a single field, 'max' or 'full', not an array of structs
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,struct('max',1,'full',{{[2 1 2],[0 1 1]}})),'single field');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,struct('mx',1)),'single field');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,struct('full',{[2 1 2],[0 1 1]})),'scalar struct');
assert_err(@() piesos_intopvar(prog,[1,1],vars,dom,[1,2],intmat12,struct('full',{{[2 1 2]}})),'with an element for each term');
end

function case_no_variables()
% No spatial variables (a constant matrix variable); variables with deps=0 are dropped
prog = new_prog({'s1','s2'});
[prog,K0] = piesos_intopvar(prog,[2,3],cell(1,0),zeros(0,2),zeros(1,0),'full',3);
assert(isequal(size(K0.params),[1,1]) && isempty(K0.vars.in),'Expected a constant functional.');
assert(isa(ker(K0,1),'dpvar') && ndv(K0,1)==6 && isequal(size(K0),[2,3]),'Wrong constant kernel.');
[~,K1] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 1],[0,2],'full',1);
assert(isequal(K1.vars.in,{'s2'}) && isequal(size(K1.params),[3,1]),'Variable with deps=0 must be dropped.');
assert(isequal(K1.dom.in,[0 1]),'Wrong domain.');
% also when intmat and degrees are given by cell
[~,K2] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[0,2],{zeros(1,0),[0,2,0,3]},{0,[1 1]});
assert(isequal(K2.vars.in,{'s2'}) && isequal(K2.dom.in,[0 2]) && numel(K2.params)==1,'Wrong functional.');
assert(ndv(K2,1)==4,'Maximal degrees over the dummy variables of s2: 2x2 monomials.');
[~,K3] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[0,2],{zeros(1,0),[0,2,0,3]},dfull({{[],[0 0;1 1;1 0;0 1]}}));
assert(isequal(K3.vars.in,{'s2'}) && numel(K3.params)==1,'Wrong functional with per-term degrees.');
assert(ndv(K3,1)==4,'The degrees of the discarded variable s1 must be dropped: s2 has 4 monomials.');
% a 1 x N cell of per-variable monomials, with an array without columns for the discarded variable
[~,K3b] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[0,2],{zeros(1,0),[0,2,0,3]},dfull({zeros(1,0),[0 0;1 1;1 0;0 1]}));
assert(isequal(K3b.vars.in,{'s2'}) && ndv(K3b,1)==4,'Monomials of s2 as 1 x N cell, without s1.');
% all variables dropped
[~,K4] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 1],[0,0],'full',2);
assert(isempty(K4.vars.in) && numel(K4.params)==1,'Expected a constant functional.');
% A single variable with deps=0 (indexing a scalar by false gives 0 x 0, not 1 x 0), for
% a generated and an explicit intmat, with maximal degrees, and with per-variable
% monomials (an element for the discarded variable is ignored)
[~,K5] = piesos_intopvar(prog,[1,1],{'s1'},[0 1],0,'full',2);
assert(isempty(K5.vars.in) && numel(K5.params)==1 && ndv(K5,1)==1,'Expected a constant functional for N=1, deps=0.');
[~,K6] = piesos_intopvar(prog,[1,1],{'s1'},[0 1],0,{zeros(1,0)},dmax({7}));
assert(isempty(K6.vars.in) && numel(K6.params)==1 && ndv(K6,1)==1,'Expected a constant functional for an explicit intmat.');
[~,K7] = piesos_intopvar(prog,[1,1],{'s1'},[0 1],0,'full',dfull({[]}));
assert(isempty(K7.vars.in) && numel(K7.params)==1 && ndv(K7,1)==1,'Expected a constant functional for full degrees.');
% the array of a discarded variable (deps=0) must have no columns, in a 1 x N cell
% of the variables (else it is not recognized) and in the 1 x N cell of a term
assert_err(@() piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[0,2],{zeros(1,0),[0,2,0,3]},dfull({5,[0 0;1 1]})),'ntrms(1) x ... x ntrms(N) cell');
assert_err(@() piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[0,2],{zeros(1,0),[0,2,0,3]},dfull({{5,[0 0;1 1]}})),'deps(i) columns');
end

function case_struct_inputs()
prog = new_prog({'s1','s2'});
[~,K1] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[1,2],'nomult',1);
[~,K2] = piesos_intopvar(prog,[1,1],struct('in',{{'s1','s2'}},'out',{cell(1,0)}),...
            struct('in',[0 1;0 2],'out',zeros(0,2)),[1,2],'nomult',1);
assert(isequal(K1.vars,K2.vars) && isequal(K1.dom,K2.dom) && isequal(K1.intmat,K2.intmat),'Struct inputs must give the same functional.');
% a single interval and a scalar deps apply to all variables
[~,K3] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 2],1,'full',0);
assert(isequal(K3.dom.in,[0 2;0 2]) && isequal(K3.depmat.in,[1 1]),'Single interval and scalar deps must apply to all variables.');
% column of variable names, and a functional with zero rows
[~,Kc] = piesos_intopvar(prog,[1,1],{'s2';'s1'},[0 1;0 2],[1,1],'full',0);
assert(isequal(Kc.vars.in,{'s1','s2'}) && isequal(Kc.dom.in,[0 2;0 1]),'Column of variables must be sorted with domains.');
[~,Kz] = piesos_intopvar(prog,[0,2],{'s1'},[0 1],2,'full',1);
assert(isequal(size(Kz),[0,2]) && isvalid(Kz),'A functional with zero rows must be valid.');
% string option
[~,K4] = piesos_intopvar(prog,[1,1],{'s1','s2'},[0 1;0 2],[1,2],"nomult",1);
assert(isequal(K4.intmat,K1.intmat),'A string option must work as a char.');
end

function case_sequential()
% Two functionals in the same program have distinct decision variables, and
% the first is unchanged by declaring the second
prog = new_prog({'s'});
[prog,K1] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',1);
n1 = prog.var.num;
[prog,K2] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',1);
assert(prog.var.num==2*n1,'Second functional must add its own decision variables.');
names1 = cellfun(@(P) P.dvarname(:),kers(K1),'UniformOutput',false);
names2 = cellfun(@(P) P.dvarname(:),kers(K2),'UniformOutput',false);
names1 = vertcat(names1{:});    names2 = vertcat(names2{:});
assert(isempty(intersect(names1,names2)),'Decision variables of different functionals must be distinct.');
end

function case_arithmetic()
% Functionals declared by piesos_intopvar (different intmats, different
% decision variables) combine: the sum applies as the sum of the functionals
rng(10);
prog = new_prog({'s'});
[prog,K1] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'lower',1);
[prog,K2] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'upper',1);
[prog,K3] = piesos_intopvar(prog,[1,1],{'s'},[0,1],2,'full',1);
d0 = rand_dvals({K1,K2,K3});
Ks = K3-2*K1+K2;
assert(isvalid(Ks),'Combined functional must be valid.');
[~,vnames] = names_of(K3);
v = rand_polyvec(vnames,1,5,2);
ev = @(K) double(apply_op(eval_intop_dpvar(K,d0),v));
assert(abs(ev(Ks)-(ev(K3)-2*ev(K1)+ev(K2)))<1e-10*max(1,abs(ev(Ks))),'Combined functional must act as the sum.');
end

function [is_ok,prog] = solve_ok(prog)
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
is_ok = (info.pinf==0 && info.dinf==0);
end

function case_sdp_recovery()
% Kop*v = F(v) for all monomials v up to degree 3 in each of s_1,s_2, with F
% the functional with kernel kappa on [a,b]^2, kappa of degree 1 in each
% variable. Needs kernels on both s_1>=s_2 and s_1<=s_2 (the same kernel):
% feasible for 'nomult' and 'full', not for 'lower' nor for constant kernels
a = 0.5;    b = 2;
kap = @(s1,s2) 1+2*s1-s2+0.5*s1.*s2;
sv = polynomial({'s_1','s_2'});
F = zeros(4,4);
for i=0:3
    for j=0:3
        F(i+1,j+1) = integral2(@(s1,s2) kap(s1,s2).*s1.^i.*s2.^j,a,b,a,b,'AbsTol',1e-13,'RelTol',1e-13);
    end
end
build = @(opt,degs) build_recovery_prog(opt,degs,[a,b],sv,F);
% (a) feasible, kernel recovered
[ok,prog,K] = build('nomult',[1 1]);
assert(ok,'The functional with kernel kappa must be feasible for ''nomult''.');
Ks = piesos_getsol(prog,K);
for t=1:2
    for pt=[0.6 0.9; 1.1 1.7; 1.9 0.7]'
        val = evalpoly(ker(Ks,t),{'s_dum1','s_dum2'},pt');
        assert(abs(val-kap(pt(1),pt(2)))<1e-5,'Kernel %d not recovered: %g vs %g.',t,val,kap(pt(1),pt(2)));
    end
end
% also on monomials not used in the constraints
v = polynomial({'s_1','s_2'});     v = v(1)^4*v(2) + 2*v(2)^2;
cnum = double(apply_op(Ks,v));
cref = integral2(@(s1,s2) kap(s1,s2).*(s1.^4.*s2+2*s2.^2),a,b,a,b,'AbsTol',1e-13,'RelTol',1e-13);
assert(abs(cnum-cref)<1e-4*max(1,abs(cref)),'Recovered functional differs on other polynomials.');
% (b) 'full' is feasible as well
ok = build('full',[1 1]);
assert(ok,'The functional must be feasible for ''full'' as well.');
% (c) infeasible with only s_1>=s_2, or with constant kernels
ok = build('lower',[1 1]);
assert(~ok,'The functional must be infeasible for ''lower''.');
ok = build('nomult',[0 0]);
assert(~ok,'The functional must be infeasible with constant kernels.');
end

function [ok,prog,K] = build_recovery_prog(opt,degs,dom,sv,F)
prog = new_prog({'s'});
[prog,K] = piesos_intopvar(prog,[1,1],{'s'},dom,2,opt,{degs});
% SeDuMi cannot solve an equality-only program: add a trivial PSD matrix
[prog,~] = sosposmatr(prog,1);
for i=0:3
    for j=0:3
        c = apply_op(K,sv(1)^i*sv(2)^j);
        prog = soseq(prog,c-F(i+1,j+1));
    end
end
[ok,prog] = solve_ok(prog);
end

function case_sdp_two_variables()
% Two spatial variables of degree 1: the functional with kernel kappa(s1,s2)
% on [a1,b1]x[a2,b2], recovered from its values on monomials s1^i s2^j
rng(9);
dom = [0 1;1 3];
kap = @(s1,s2) 2+s1.*s2+s2.^2;
prog = new_prog({'s1','s2'});
[prog,K] = piesos_intopvar(prog,[1,1],{'s1','s2'},dom,[1,1],'full',{1,2});
[prog,~] = sosposmatr(prog,1);
sv = polynomial({'s1','s2'});
for i=0:2
    for j=0:3
        Fij = integral2(@(s1,s2) kap(s1,s2).*s1.^i.*s2.^j,dom(1,1),dom(1,2),dom(2,1),dom(2,2),'AbsTol',1e-13,'RelTol',1e-13);
        prog = soseq(prog,apply_op(K,sv(1)^i*sv(2)^j)-Fij);
    end
end
[ok,prog] = solve_ok(prog);
assert(ok,'The functional must be feasible.');
Ks = piesos_getsol(prog,K);
for pt=[0.2 1.5; 0.7 2.2; 0.9 2.9]'
    val = evalpoly(ker(Ks,1),{'s1_dum1','s2_dum1'},pt');
    assert(abs(val-kap(pt(1),pt(2)))<1e-5,'Kernel not recovered: %g vs %g.',val,kap(pt(1),pt(2)));
end
% With kernels of degree 0 in s2, infeasible
prog = new_prog({'s1','s2'});
[prog,K] = piesos_intopvar(prog,[1,1],{'s1','s2'},dom,[1,1],'full',{1,0});
[prog,~] = sosposmatr(prog,1);
for i=0:2
    for j=0:3
        Fij = integral2(@(s1,s2) kap(s1,s2).*s1.^i.*s2.^j,dom(1,1),dom(1,2),dom(2,1),dom(2,2),'AbsTol',1e-13,'RelTol',1e-13);
        prog = soseq(prog,apply_op(K,sv(1)^i*sv(2)^j)-Fij);
    end
end
assert(~solve_ok(prog),'Kernels of degree 0 in s2 cannot represent kappa.');
end

function case_errors()
prog = new_prog({'s1','s2'});
v = {'s1','s2'};    D = [0 1;0 1];
assert_err(@() piesos_intopvar(prog,[1,2,3],v,D,1,'full',0),'Matrix dimensions');
assert_err(@() piesos_intopvar(prog,[1,-2],v,D,1,'full',0),'Matrix dimensions');
assert_err(@() piesos_intopvar(prog,[1,1],{'s1','s1'},D,1,'full',0),'distinct');
assert_err(@() piesos_intopvar(prog,[1,1],'s1',D,1,'full',0),'Variables must be specified');
assert_err(@() piesos_intopvar(prog,[1,1],struct('out',{{}}),D,1,'full',0),'No input variables');
assert_err(@() piesos_intopvar(prog,[1,1],struct('in',{v},'out',{v}),D,1,'full',0),'output variables');
assert_err(@() piesos_intopvar(prog,[1,1],v,[0 1 2;0 1 2],1,'full',0),'upper and lower boundary');
assert_err(@() piesos_intopvar(prog,[1,1],v,[0 1;0 1;0 1],1,'full',0),'Number of domains');
assert_err(@() piesos_intopvar(prog,[1,1],v,[1 0;0 1],1,'full',0),'Lower boundary');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1 -1],'full',0),'nonnegative integers');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1 1.5],'full',0),'nonnegative integers');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1 1 1],'full',0),'1 x N array');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'everything',0),'Option for');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,3,0),'cell or as char');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,{[0,2]},0),'1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[2,1],{[0,2],[0,2]},0),'2*q columns');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,{zeros(0,2),[0,2]},0),'At least one interval');
% an invalid intmat (variables integrated over themselves) is rejected before
% any decision variable is declared
assert_err(@() piesos_intopvar(prog,[1,1],{'s1'},[0 1],2,{[0,1,0,1]},0),'not properly specified');
% degrees
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'full',[1 2]),'scalar, or as 1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'full',{1}),'1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'full',{1,-1}),'nonnegative integers');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'full',{1,1.5}),'nonnegative integers');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'full','a'),'scalar, or as 1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[2,1],'full',{[1 2 3],0}),'scalar or 1 x d array');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[2,1],'full',{ones(2,3),0}),'Use the field ''full''');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[2,1],'full',{[0 1;1 1],0}),'Use the field ''full''');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,1,'full',{{0,0}}),'1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1,1],{[0,2],[0,2]},dfull({{0,0},{0}})),'with an element for each term');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1,1],{[0,2],[0,2]},dfull({{0}})),'1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1,1],{[0,2],[0,2]},[0 0]),'scalar, or as 1 x N cell');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1,1],{[0,2],[0,2]},dfull(0)),'Full degrees must be specified');
assert_err(@() piesos_intopvar(prog,[1,1],v,D,[1,1],{[0,2],[0,2]},dmax([0 0])),'scalar, or as 1 x N cell');
% a 1 x N cell for a cell intmat is the maximal degrees of the variables
[~,Km] = piesos_intopvar(prog,[1,1],v,D,[1,1],{[0,2],[0,2]},{1,1});
assert(ndv(Km,1)==4,'Maximal degrees (1,1): 4 monomials.');
% a valid call after errors
[~,K] = piesos_intopvar(prog,[1,1],v,D,[1,1],'full',0);
assert(isvalid(K),'Valid call must work after errors.');
end

function Kk = kers(K)
% Kernels of K (each 'dpvar' kernel in the decision variables it uses)
Kk = kernels(K);
end

function P = ker(K,k)
% Kernel k of K, in the decision variables that it uses
Kk = kernels(K);
P = Kk{k};
end

function n = ndv(K,k)
% Number of decision variables used by kernel k of K
P = ker(K,k);
if isa(P,'dpvar')
    n = numel(P.dvarname);
else
    n = 0;
end
end
