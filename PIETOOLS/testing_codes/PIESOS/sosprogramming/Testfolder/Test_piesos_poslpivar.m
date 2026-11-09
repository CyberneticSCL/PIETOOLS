function Test_piesos_poslpivar()
% TEST_PIESOS_POSLPIVAR tests, using assertions, piesos_poslpivar,
%   [prog,Pmat,Zop] = ...
%               piesos_poslpivar(prog,dims,vars,dom,depmat,degs,opts),
% which declares Pop = Zop'*Pmat*Zop, with Zop a basis of tensor-PI
% operators whose kernels depend only on the dummy variables, and Pmat a
% q x q cell of SOS matrices in the primary variables (one row/column per
% row of Zop.ops).
%
% Ground truths are independent of the function under test:
%  - Zop: the number of rows is counted by hand from the monomial bases. A
%    3-PI direction with dummy degree r gives 1 multiplier row (no dummy
%    monomial) and (r+1) rows for each of the lower and upper integral;
%    the counts multiply over the directions. Every factor must have ZL=0
%    (no primary monomial), and ZR the dummy monomials (0:r)';
%  - Pmat: for each psatz term, a random positive semidefinite Gram matrix
%    Qk is assigned to the decision variables, through the names returned
%    by a direct call of sosquadvar on a fresh program with monomial
%    vectors Z built by hand. Pmat(s), evaluated at random points s, is
%    compared with
%        sum_k g_k(s) * (I x Z(s))' * Qk * (I x Z(s)),
%    with Z and the multipliers g_k written out by hand in each case. The
%    declared decision variables are compared with those of the direct
%    call;
%  - degrees specified in any of the formats of monopmat are compared with
%    the interpretation by monopmat itself (an unsplit call, whose ZL
%    holds the monomials in the output variables);
%  - SDP: the declared Pmat is used in actual programs (sossolve, SeDuMi),
%    where Pmat = target is feasible only for targets that are SOS
%    matrices (or, with a psatz multiplier, positive on the domain), and
%    jointly PSD over all blocks;
%  - end to end: <Zop*x, Pmat*Zop*x>_{L2} is computed by polynomial
%    arithmetic (apply_op, then integration) for arbitrary polynomial x,
%    and compared with the sum of squares
%    sum_k int g_k*|Rk*(I x Z)*Zop*x|^2, with Qk = Rk'*Rk, built from the
%    Cholesky factors without using Pmat. This also shows that
%    Pop = Zop'*Pmat*Zop is positive.
% The operator Pop itself is not returned (see piesos_poslpivar).
%
% States depending on different variables (mixed dependence) are all
% mapped to the functions on the common domain of the variables on which
% some state depends (see piesos_poslpivar), so that Pop is positive for
% the L2 inner product on that domain. The mixed cases (case_mixed_*, and
% the extensions of case_degs_formats, case_positive_form) check that
%  - the row counts and monomials of Zop are those of each state's own
%    input variables, with ZL=0 and output variables covering all used
%    variables;
%  - Pmat blocks and psatz multipliers are functions of ALL used
%    variables, also for a block whose state does not depend on all of
%    them;
%  - <Zop*x,Pmat*Zop*x> is the same by polynomial arithmetic (apply_op),
%    by the sum of squares of the Cholesky factors, and through the route
%    used in the programs, innerprod(Zop*Zx,Pmat), with Zx the 'polyopvar'
%    states;
%  - in an SDP, a psatz multiplier in a variable on which the first state
%    does not depend still acts on its block (target g(s2)*I on all
%    blocks).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_poslpivar
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
% DJ (CC), 10/21/2026: 3 x 1 'exclude' is OR over the variables of a state
%                        (nrows_nd, case_mixed_exclude_sep);

% NOTES
% - The SDP cases (case_sdp_feasibility,
%   case_mixed_sdp) call sossolve with SeDuMi;
% - The innerprod route (innerprod_functional) is used for scalar states
%   in <=2 variables only: innerprod fails for mixed states in 3
%   variables;
% - case_unsorted_variables uses dummy degree 0 and does not check Zop:
%   for a positive degree, monopmat returns a different operator than for
%   alphabetical order (warning 'sopvar:noncanonicalMultiplier');
% - nrows_nd: monopmat applies a 3-vector exclude along the variables of
%   each state, and excludes a term if its kind is excluded in any of them.

tests = {@case_split_1D, @case_defaults, @case_degs_formats, ...
         @case_struct_inputs, @case_2D, @case_3D, @case_two_states, ...
         @case_distributed_deg2, @case_exclude_sep, @case_psatz, ...
         @case_psatz_logical, @case_unused_variable, ...
         @case_unsorted_variables, @case_sequential, @case_no_variables, ...
         @case_positive_form, @case_sdp_feasibility, @case_errors, ...
         @case_mixed_structure, @case_mixed_exclude_sep, ...
         @case_mixed_psatz, @case_mixed_sdp};
nfail = 0;
for i=1:numel(tests)
    tstart = tic;
    try
        tests{i}();
        fprintf('PASS  %s (%.1f s)\n',func2str(tests{i}),toc(tstart));
    catch err
        nfail = nfail+1;
        fprintf('FAIL  %s\n      %s\n',func2str(tests{i}),err.message);
        for k=1:min(3,numel(err.stack))
            fprintf('        in %s, line %d\n',err.stack(k).name,err.stack(k).line);
        end
    end
end
if nfail>0
    error('Test_piesos_poslpivar: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_poslpivar: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function prog = new_prog(varnames)
% Fresh program in the given spatial variables
N = numel(varnames);
x1 = polyopvar('x1',polynomial(varnames(:)),[zeros(N,1),(1:N)']);
prog = piesos_program(x1);
end

function m = nrows_3pi(r,excl)
% Rows of the basis in one 3-PI direction, for dummy degree r (kernel
% monomials (0:r)', no primary monomial): excl = logical [mult,lower,upper]
% of excluded terms. The multiplier has no dummy monomial, so a single row.
m = ~excl(1)*1 + ~excl(2)*(r+1) + ~excl(3)*(r+1);
end

function z = zvec_num(E,s)
% Numeric monomial vector Z(s) from the exponents E{k} of the variables k,
% with the first variable slowest (Kronecker order)
z = 1;
for k=1:numel(E)
    z = kron(z,s(k).^E{k}(:));
end
end

function Z = zvec_poly(E,varnames)
% Polynomial monomial vector with the same order as zvec_num, built from
% the mixed-radix expansion of the index, independently of the code
% under test
if isempty(E)
    Z = polynomial(1);
    return
end
nE = cellfun(@numel,E);
nZ = prod(nE);
pv = polynomial(varnames(:));
Z = cell(nZ,1);
for idx=0:nZ-1
    r = idx;
    mon = polynomial(1);
    for k=numel(E):-1:1
        ik = mod(r,nE(k))+1;
        r = floor(r/nE(k));
        mon = mon*pv(k)^E{k}(ik);
    end
    Z{idx+1} = mon;
end
Z = vertcat(Z{:});
end

function Q = random_psd(n)
R = randn(n);
Q = R'*R/n + 0.1*eye(n);
Q = (Q+Q')/2;
end

function v = evalpoly(P,varnames,vals)
% Evaluate the polynomial matrix P at the point vals of the variables
% varnames, using the fields of the 'polynomial' object directly
if isa(P,'double')
    v = P;
    return
end
C = full(P.coefficient);
deg = full(P.degmat);
mon = ones(size(deg,1),1);
for k=1:numel(P.varname)
    idx = find(strcmp(varnames,P.varname{k}));
    assert(~isempty(idx),'Pmat depends on unexpected variable %s.',P.varname{k});
    mon = mon.*vals(idx).^deg(:,k);
end
v = reshape(C.'*mon,P.matdim);
end

function E = ref_exponents(dims,varnames,dom,depmat_in,degs,opts)
% Exponents of the monomials in the primary variables, as interpreted by
% monopmat itself: ZL of an unsplit call with every state mapped to the
% variables on which some state depends. E{i}{k}: kth USED variable
vars = struct('in',{varnames},'out',{varnames});
doms = struct('in',dom,'out',dom);
q = numel(dims);
used = any(depmat_in>0,1);
dep = struct('in',depmat_in,'out',repmat(used,q,1));
if nargin<6, opts = struct(); end
Zref = monopmat(dims,vars,doms,dep,degs,opts);
vused = varnames(used);
E = cell(1,q);
for i=1:q
    sop = Zref.ops{i,i}.ops{1};
    E{i} = cell(1,numel(vused));
    for k=1:numel(vused)
        E{i}{k} = sop.ZL{strcmp(sop.vars.out,vused{k})};
    end
end
end

function e = exps_by_name(vals,names,varnames)
% Reorder the cell vals, indexed as the cellstr names, as the cellstr varnames
% ([] for a variable that is not in names)
e = cell(1,numel(varnames));
for k=1:numel(varnames)
    idx = find(strcmp(names,varnames{k}));
    if ~isempty(idx), e{k} = vals{idx}; end
end
end

function check_Zop(Zop,m_exp,ZR_exp,varnames,vout)
% The kernels have no primary monomial (ZL=0) and the dummy monomials ZR,
% for each block row: ZR_exp{i} = 1 x N cell of exponents (first factor),
% [] for a variable on which the state does not depend. If vout is given, a
% block with a single factor must have exactly these output variables (all
% variables on which some state depends)
q = numel(m_exp);
assert(isequal(size(Zop.ops),[q,q]),'Zop must have q x q blocks.');
assert(isequal(Zop.dims.out(:),m_exp(:)),'Row counts of Zop blocks: expected %s, got %s.',...
    mat2str(m_exp(:)'),mat2str(Zop.dims.out(:)'));
assert(size(Zop,1)==sum(m_exp),'Total row count of Zop is wrong.');
for i=1:q
    facs = Zop.ops{i,i}.ops;
    for f=1:numel(facs)
        sop = facs{f};
        assert(all(cellfun(@(z) isequal(z(:),0),sop.ZL)),...
            'Block %d, factor %d: kernel depends on primary variables (ZL~=0).',i,f);
    end
    sop = Zop.ops{i,i}.ops{1};
    ZR = exps_by_name(sop.ZR,sop.vars.in,varnames);
    for k=1:numel(varnames)
        assert(isequal(ZR{k}(:),ZR_exp{i}{k}(:)),'Block %d: wrong dummy monomials in variable %s.',i,varnames{k});
    end
    if nargin>=5 && numel(facs)==1
        assert(isempty(setxor(sop.vars.out(:),vout(:))),...
            'Block %d: output variables must be %s.',i,strjoin(vout,','));
    end
end
end

function check_sos(prog,Pmat,Zop,spec,prog_fresh)
% Pmat = sum_k g_k*(I x Z)'*Qk*(I x Z) for random PSD Gram matrices Qk
% spec fields: varnames (1 x N), dom (N x 2), E (1 x q cell of 1 x N cell of
% exponents of Z_i), g (1 x T cell of handles g(s), in the order of declaration)
q = numel(spec.E);
T = numel(spec.g);
N = numel(spec.varnames);
assert(iscell(Pmat) && isequal(size(Pmat),[q,q]),'Pmat must be a q x q cell.');
m = zeros(q,1);     nz = zeros(q,1);
Zs = cell(q,1);
for i=1:q
    m(i) = size(Zop.ops{i,i},1);
    nz(i) = prod(cellfun(@numel,spec.E{i}));
    Zs{i} = zvec_poly(spec.E{i},spec.varnames);
    for j=1:q
        assert(isequal(size(Pmat{i,j}),[size(Zop.ops{i,i},1),size(Zop.ops{j,j},1)]),...
            'Pmat{%d,%d} has wrong size.',i,j);
        assert(isa(Pmat{i,j},'dpvar'),'Pmat{%d,%d} must be a dpvar.',i,j);
    end
end
% Direct declarations in a fresh program, in the order of the terms
prog0 = prog_fresh;
Qn = cell(1,T);
for t=1:T
    [prog0,~,Qn{t}] = sosquadvar(prog0,Zs,Zs,m,m,'pos');
end
assert(prog.var.num==prog_fresh.var.num+T && isequal(prog.var.idx,prog0.var.idx),...
    'Declared decision variables differ from %d direct sosquadvar call(s).',T);
assert(all(strcmp(prog.var.type(prog_fresh.var.num+1:end),'sos')),'All new variables must be declared as sos.');
% Assign random PSD Gram matrices
nQ = sum(m.*nz);
off = [0;cumsum(m.*nz)];
names = cell(0,1);      vals = zeros(0,1);
Qnum = cell(1,T);
for t=1:T
    Qnum{t} = random_psd(nQ);
    for i=1:q
        for j=1:q
            nm = Qn{t}{i,j};
            vl = Qnum{t}(off(i)+1:off(i+1),off(j)+1:off(j+1));
            names = [names;nm(:)];          %#ok<AGROW>
            vals = [vals;vl(:)];            %#ok<AGROW>
        end
    end
end
[names,iu] = unique(names,'stable');
vals = vals(iu);
d0map = containers.Map(names,num2cell(vals));
% Gram matrices assigned: polynomial blocks, evaluated once (eval_dpvar_at is
% slow for large Gram matrices)
Pevals = cell(q,q);
for i=1:q
    for j=1:q
        Pevals{i,j} = eval_dpvar_at(Pmat{i,j},d0map);
    end
end
% Compare at random points in the domain
for p=1:3
    s0 = spec.dom(:,1)' + rand(1,N).*(spec.dom(:,2)-spec.dom(:,1))';
    Pfull = zeros(sum(m));
    for i=1:q
        for j=1:q
            Pe = Pevals{i,j};
            Pnum = evalpoly(Pe,spec.varnames,s0);
            Pexp = zeros(m(i),m(j));
            zi = zvec_num(spec.E{i},s0);        zj = zvec_num(spec.E{j},s0);
            for t=1:T
                Qij = Qnum{t}(off(i)+1:off(i+1),off(j)+1:off(j+1));
                Pexp = Pexp + spec.g{t}(s0)*kron(eye(m(i)),zi)'*Qij*kron(eye(m(j)),zj);
            end
            assert(norm(Pnum-Pexp,'fro')<=1e-8*(1+norm(Pexp,'fro')),...
                'Pmat{%d,%d}(s) differs from the SOS expression (error %g).',i,j,norm(Pnum-Pexp,'fro'));
            r1 = sum(m(1:i-1))+1;   c1 = sum(m(1:j-1))+1;
            Pfull(r1:r1+m(i)-1,c1:c1+m(j)-1) = Pnum;
        end
    end
    assert(norm(Pfull-Pfull','fro')<=1e-8*(1+norm(Pfull,'fro')),'Pmat(s) is not symmetric.');
    assert(min(eig((Pfull+Pfull')/2))>=-1e-8*(1+norm(Pfull,'fro')),...
        'Pmat(s) is not positive semidefinite at a point of the domain.');
end
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

function E = uniform_E(q,N,degs)
% Same exponents (0:degs(k))' in variable k, for each of q blocks
E = cell(1,q);
for i=1:q
    E{i} = cell(1,N);
    for k=1:N
        E{i}{k} = (0:degs(k))';
    end
end
end

function g = no_multiplier()
g = {@(s) 1};
end


%% ---- cases ------------------------------------------------------------
function case_split_1D()
% Single scalar state, one variable: dummy degree r (kernels of Zop) and
% primary degree d (Pmat of degree 2d) are set independently
rng(1);
for r=[0,1,3]
    for d=[0,1,2]
        prog0 = new_prog({'s1'});
        [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[0.5 2],1,struct('in',r,'out',d),struct());
        m = nrows_3pi(r,false(1,3));
        check_Zop(Zop,m,{{(0:r)'}},{'s1'});
        spec = struct('varnames',{{'s1'}},'dom',[0.5 2],'E',{uniform_E(1,1,d)},'g',{no_multiplier()});
        check_sos(prog,Pmat,Zop,spec,prog0);
        % Pmat has degree 2d in the primary variable, and size m x m
        assert(full(max([0;Pmat{1}.degmat(:)]))==2*d,'Degree of Pmat should be 2d=%d.',2*d);
        % decision variables: one Gram matrix of size m*(d+1)
        assert(prog.var.idx{2}-prog.var.idx{1}==(m*(d+1))^2,'Wrong number of decision variables.');
    end
end
end

function case_defaults()
% Omitted degs and opts: degree 0 in the primary and dummy variables
rng(2);
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[0 1],1);
check_Zop(Zop,3,{{0}},{'s1'});
spec = struct('varnames',{{'s1'}},'dom',[0 1],'E',{uniform_E(1,1,0)},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(isempty(Pmat{1}.degmat) || all(Pmat{1}.degmat(:)==0),'Pmat must be constant for degree 0.');
% Scalar degs: the same degree for the primary and dummy variables
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[0 1],1,2);
check_Zop(Zop,nrows_3pi(2,false(1,3)),{{(0:2)'}},{'s1'});
spec.E = uniform_E(1,1,2);
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_degs_formats()
% All formats of degs accepted by monopmat: the primary monomials Z(s) must
% be the ones monopmat itself assigns to the output variables (ZL of an
% unsplit call), and the dummy monomials ZR those of the input
rng(3);
dom1 = [0.5 2];                 dom2 = [0.5 2; -1 0.7];
cases = {};
% 1 variable
cases{end+1} = {1,{'s1'},dom1,1,2};                                       % scalar
cases{end+1} = {1,{'s1'},dom1,1,struct('in',1,'out',2)};                  % struct
cases{end+1} = {1,{'s1'},dom1,1,struct('in',{{1}},'out',{{2}})};          % per-state cells
cases{end+1} = {1,{'s1'},dom1,1,struct('in',{{{[0;1;2]}}},'out',{{{[0;2;3]}}})};  % monomial lists
cases{end+1} = {1,{'s1'},dom1,1,{{[0;1]}}};                               % list for in and out
cases{end+1} = {1,{'s1'},dom1,1,struct('in',2,'out',{{{[0;1;4]}}})};      % array in, list out
% 2 variables
cases{end+1} = {1,{'s1','s2'},dom2,[1 1],[2 1]};                          % per-variable degrees
cases{end+1} = {1,{'s1','s2'},dom2,[1 1],struct('in',[1 2],'out',[2 0])};
cases{end+1} = {1,{'s1','s2'},dom2,[1 1],struct('in',1,'out',[0 2])};
cases{end+1} = {1,{'s1','s2'},dom2,[1 1],struct('in',{{{[0;1],[0;2]}}},'out',{{{[0;3],[0;1;2]}}})};
% 2 states, own degrees
cases{end+1} = {[1;1],{'s1'},dom1,[1;1],struct('in',{{1,2}},'out',{{2,1}})};
cases{end+1} = {[1;1],{'s1'},dom1,[1;1],struct('in',{{1,2}},'out',{{{[0;2]},{[0;1;3]}}})};
% 2 states of different dependence (mixed): all formats, output monomials in
% both variables for both states (the output space is common), input monomials
% only in the variables a state depends on
dmix = [1 0;1 1];
cases{end+1} = {[1;1],{'s1','s2'},dom2,dmix,1};
cases{end+1} = {[1;1],{'s1','s2'},dom2,dmix,struct('in',[0 1],'out',[2 1])};
cases{end+1} = {[1;1],{'s1','s2'},dom2,dmix,struct('in',{{0,[1 0]}},'out',{{2,[1 2]}})};
cases{end+1} = {[1;1],{'s1','s2'},dom2,dmix,...
    struct('in',{{{[0;1]},{[0;1],[0;2]}}},'out',{{{[0;2],[0;1]},{[0;1;3],[0;3]}}})};
cases{end+1} = {[1;2],{'s1','s2'},dom2,[1 0;0 1],struct('in',{{1,2}},'out',{{[2 0],[0 1]}})};
cases{end+1} = {[1;1],{'s1','s2'},dom2,[0 0;1 1],struct('in',[1 1],'out',{{1,[2 1]}})};
cases{end+1} = {[1;1;1],{'s1','s2'},dom2,[1 0;0 0;0 1],struct('in',[1 1],'out',{{1,2,[1 3]}})};
for c=1:numel(cases)
    [dims,vn,dm,depmat,degs] = cases{c}{:};
    q = numel(dims);
    used = any(depmat>0,1);
    prog0 = new_prog(vn);
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,dims,vn,dm,depmat,degs,struct());
    E = ref_exponents(dims,vn,dm,depmat,degs);
    % dummy monomials and row counts, from the unsplit reference
    if isstruct(degs), din = degs.in; else, din = degs; end
    Zref = monopmat(dims,struct('in',{vn},'out',{vn}),struct('in',dm,'out',dm),...
        struct('in',depmat,'out',repmat(used,q,1)),struct('in',{din},'out',0),struct());
    for i=1:q
        facs = Zop.ops{i,i}.ops;     facr = Zref.ops{i,i}.ops;
        assert(isequal(size(Zop.ops{i,i}),size(Zref.ops{i,i})),'Case %d: block %d has wrong size.',c,i);
        for f=1:numel(facs)
            assert(isequal(facs{f}.ZL,facr{f}.ZL) && isequal(facs{f}.ZR,facr{f}.ZR),...
                'Case %d: kernel monomials differ from monopmat with degs.out=0.',c);
        end
        assert(all(cellfun(@(z) isequal(z(:),0),facs{1}.ZL)),'Case %d: ZL must be 0.',c);
    end
    spec = struct('varnames',{vn(used)},'dom',dm(used,:),'E',{E},'g',{no_multiplier()});
    check_sos(prog,Pmat,Zop,spec,prog0);
end
end

function case_struct_inputs()
% struct-valued dims/vars/dom/depmat accepted; vector-valued state (n=2)
rng(4);
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,struct('in',2,'out',2),...
    struct('in',{{'s1'}}),struct('in',[0 1]),struct('in',1),struct('in',1,'out',1),struct());
m = 2*nrows_3pi(1,false(1,3));
check_Zop(Zop,m,{{(0:1)'}},{'s1'});
spec = struct('varnames',{{'s1'}},'dom',[0 1],'E',{uniform_E(1,1,1)},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
% Equivalent to the direct specification
prog1 = new_prog({'s1'});
[prog1,Pmat1] = piesos_poslpivar(prog1,2,{'s1'},[0 1],1,1);
assert(isequal(prog1.var.idx,prog.var.idx) && isequal(Pmat1{1}.C,Pmat{1}.C),'struct input differs from direct input.');
% Column cellstr of variables
prog1 = new_prog({'s1'});
[prog1,Pmat1] = piesos_poslpivar(prog1,2,{'s1'}',[0 1],1,1);
assert(isequal(prog1.var.idx,prog.var.idx) && isequal(Pmat1{1}.C,Pmat{1}.C),'column variable list differs from row list.');
end

function case_2D()
% Two spatial variables, mixed degrees in the dummy and primary variables
rng(5);
dom = [0.5 2; -1 0.7];
dcases = {{[1 1],[1 1]},{[2 1],[1 0]},{[1 2],[0 1]}};
for c=1:numel(dcases)
    r = dcases{c}{1};   d = dcases{c}{2};
    prog0 = new_prog({'s1','s2'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1','s2'},dom,[1 1],struct('in',r,'out',d),struct());
    m = nrows_3pi(r(1),false(1,3))*nrows_3pi(r(2),false(1,3));
    check_Zop(Zop,m,{{(0:r(1))',(0:r(2))'}},{'s1','s2'});
    spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(1,2,d)},'g',{no_multiplier()});
    check_sos(prog,Pmat,Zop,spec,prog0);
    assert(full(max(Pmat{1}.degmat(:,strcmp(Pmat{1}.varname,'s1'))))==2*d(1),'Degree in s1 should be 2*d(1).');
    assert(full(max(Pmat{1}.degmat(:,strcmp(Pmat{1}.varname,'s2'))))==2*d(2),'Degree in s2 should be 2*d(2).');
end
end

function case_3D()
% Three spatial variables, degree 0 in the primary variables (so that the
% number of decision variables remains small), and one primary degree
rng(6);
dom = [0.5 2; -1 0.7; 1 3];
for dc={[0 0 0],[1 0 0]}
    d = dc{1};
    prog0 = new_prog({'s1','s2','s3'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1','s2','s3'},dom,[1 1 1],struct('in',[1 1 0],'out',d),struct());
    m = nrows_3pi(1,false(1,3))^2*nrows_3pi(0,false(1,3));
    check_Zop(Zop,m,{{(0:1)',(0:1)',0}},{'s1','s2','s3'});
    spec = struct('varnames',{{'s1','s2','s3'}},'dom',dom,'E',{uniform_E(1,3,d)},'g',{no_multiplier()});
    check_sos(prog,Pmat,Zop,spec,prog0);
end
end

function case_two_states()
% Two states with their own degrees; the second is vector-valued (n=2).
% Pmat{i,j} couples the ith and jth block row of Zop.ops
rng(7);
dom = [0.5 2];
prog0 = new_prog({'s1'});
degs = struct('in',{{1,2}},'out',{{2,1}});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;2],{'s1'},dom,[1;1],degs,struct());
m = [nrows_3pi(1,false(1,3)); 2*nrows_3pi(2,false(1,3))];
check_Zop(Zop,m,{{(0:1)'},{(0:2)'}},{'s1'});
spec = struct('varnames',{{'s1'}},'dom',dom,'E',{{{(0:2)'},{(0:1)'}}},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
% Pmat is symmetric as a block matrix, Pmat{1,2} = Pmat{2,1}'
for i=1:2
    for j=1:2
        assert(isequal(size(Pmat{i,j}),[m(i),m(j)]),'Pmat{%d,%d} has wrong size.',i,j);
    end
end
D = combine(Pmat{1,2}-Pmat{2,1}.');
assert(isempty(D.C) || all(abs(D.C(:))==0),'Pmat{1,2} must equal Pmat{2,1}''.');
% A single Gram matrix couples all blocks: exactly one declared variable
assert(prog.var.num==1,'Expected a single declared SOS variable.');
% Three states in two variables
prog0 = new_prog({'s1','s2'});
dom2 = [0.5 2; -1 0.7];
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1;1],{'s1','s2'},dom2,[1 1;1 1;1 1],struct('in',0,'out',[1 0]),struct());
assert(isequal(size(Pmat),[3,3]) && isequal(size(Zop.ops),[3,3]),'Pmat and Zop.ops must be 3 x 3.');
spec = struct('varnames',{{'s1','s2'}},'dom',dom2,'E',{uniform_E(3,2,[1 0])},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_distributed_deg2()
% depmat=2: product of two operator factors, the cumulative dummy degree r
% split over the factors as [ceil(r/2),r-ceil(r/2)]; primary monomials of
% degree d on the first factor (so one set of monomials in Pmat), and each
% factor has (1+2(rk+1)) rows, so m = prod_k(1+2(rk+1))
rng(8);
dom = [0.5 2];
for r=[1,2,3]
    rk = [ceil(r/2),r-ceil(r/2)];
    d = 1;
    prog0 = new_prog({'s1'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},dom,2,struct('in',r,'out',d),struct());
    m = prod(1+2*(rk+1));
    assert(size(Zop,1)==m,'depmat=2, r=%d: m=%d expected, got %d.',r,m,size(Zop,1));
    facs = Zop.ops{1,1}.ops;
    assert(numel(facs)==2,'Expected two factors.');
    for f=1:2
        assert(all(cellfun(@(z) isequal(z(:),0),facs{f}.ZL)),'ZL must be 0 in factor %d.',f);
        assert(isequal(facs{f}.ZR{1}(:),(0:rk(f))'),'Wrong dummy monomials in factor %d.',f);
    end
    spec = struct('varnames',{{'s1'}},'dom',dom,'E',{uniform_E(1,1,d)},'g',{no_multiplier()});
    check_sos(prog,Pmat,Zop,spec,prog0);
end
end

function case_exclude_sep()
% Excluded multiplier (as used for W in the examples), excluded integrals,
% and separability
rng(9);
dom = [0.5 2];
excls = {[1 0 0],[0 1 0],[0 0 1],[0 1 1]};
for k=1:numel(excls)
    for r=[1,2]
        prog0 = new_prog({'s1'});
        opts = struct('exclude',logical(excls{k}(:)));
        [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},dom,1,struct('in',r,'out',1),opts);
        m = nrows_3pi(r,logical(excls{k}));
        assert(size(Zop,1)==m,'exclude=%s, r=%d: m=%d expected, got %d.',mat2str(excls{k}),r,m,size(Zop,1));
        spec = struct('varnames',{{'s1'}},'dom',dom,'E',{uniform_E(1,1,1)},'g',{no_multiplier()});
        check_sos(prog,Pmat,Zop,spec,prog0);
    end
end
% Separable: the upper integral shares the parameter of the lower one
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},dom,1,struct('in',2,'out',1),struct('sep',true));
assert(size(Zop,1)==nrows_3pi(2,[false false true]),'sep: wrong row count.');
spec = struct('varnames',{{'s1'}},'dom',dom,'E',{uniform_E(1,1,1)},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_psatz()
% psatz multipliers gj = (s-a)^j*(b-s)^(q-j), j=0..q in each variable; all
% combinations; terms ordered by linear index, with the first variable
% fastest. The monomials in si are cut in degree by ceil(qi/2).
rng(10);
a = 0.5;    b = 2;
% psatz=0 is the same as no psatz
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',1,'out',2),struct('psatz',0));
spec = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,2)},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
% 1D, order q=1,2,3: terms j=0..q, primary degree d reduced by ceil(q/2)
for q=1:3
    d = 4;
    prog0 = new_prog({'s1'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',1,'out',d),struct('psatz',q));
    g = cell(1,q+1);
    for j=0:q
        g{j+1} = @(s) (s(1)-a)^j*(b-s(1))^(q-j);
    end
    dred = d-ceil(q/2);
    spec = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,dred)},'g',{g});
    check_sos(prog,Pmat,Zop,spec,prog0);
    assert(full(max(Pmat{1}.degmat(:)))<=2*d,'Degree of Pmat exceeds 2*d.');
end
% 2D, orders [1 2]: 2*3 terms, first variable fastest, vector psatz
dom = [a b; -1 0.7];
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1','s2'},dom,[1 1],struct('in',[1 0],'out',[2 3]),struct('psatz',[1 2]));
g = cell(1,6);
for l=1:6
    [j1,j2] = ind2sub([2,3],l);     j1 = j1-1;      j2 = j2-1;
    g{l} = @(s) (s(1)-dom(1,1))^j1*(dom(1,2)-s(1))^(1-j1)*(s(2)-dom(2,1))^j2*(dom(2,2)-s(2))^(2-j2);
end
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(1,2,[2-1,3-1])},'g',{g});
check_sos(prog,Pmat,Zop,spec,prog0);
% scalar psatz: the same order in all variables
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1','s2'},dom,[1 1],struct('in',[0 0],'out',[1 1]),struct('psatz',1));
g = cell(1,4);
for l=1:4
    [j1,j2] = ind2sub([2,2],l);     j1 = j1-1;      j2 = j2-1;
    g{l} = @(s) (s(1)-dom(1,1))^j1*(dom(1,2)-s(1))^(1-j1)*(s(2)-dom(2,1))^j2*(dom(2,2)-s(2))^(1-j2);
end
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(1,2,[0 0])},'g',{g});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_psatz_logical()
% Logical psatz array: only the selected multipliers, array size qi+1 along
% variable i; row vector accepted for one variable
rng(11);
a = 0.5;    b = 2;
% 1D, q=2 (size 3), terms j=0 and j=2: (b-s)^2 and (s-a)^2
for psatz = {[true;false;true],[true false true]}
    prog0 = new_prog({'s1'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',0,'out',3),struct('psatz',psatz{1}));
    g = {@(s) (b-s(1))^2, @(s) (s(1)-a)^2};
    spec = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,3-1)},'g',{g});
    check_sos(prog,Pmat,Zop,spec,prog0);
end
% 2D, 2x3 array: terms by linear index (first variable fastest)
dom = [a b; -1 0.7];
L = false(2,3);     L(1,1) = true;      L(2,2) = true;      L(1,3) = true;
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1','s2'},dom,[1 1],struct('in',[0 0],'out',[2 2]),struct('psatz',L));
lin = find(L)';
g = cell(1,numel(lin));
for t=1:numel(lin)
    [j1,j2] = ind2sub([2,3],lin(t));    j1 = j1-1;      j2 = j2-1;
    g{t} = @(s) (s(1)-dom(1,1))^j1*(dom(1,2)-s(1))^(1-j1)*(s(2)-dom(2,1))^j2*(dom(2,2)-s(2))^(2-j2);
end
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(1,2,[2-1,2-1])},'g',{g});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(prog.var.num==3,'Expected one SOS variable for each selected multiplier.');
% Logical array with only the first element true: no multiplier
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,1,struct('psatz',true));
spec = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,1)},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_unused_variable()
% A variable on which no state depends does not appear in Zop or Pmat, and a
% scalar psatz is applied to the other variable only
rng(12);
a = 0.5;    b = 2;
dom = [a b; -1 0.7];
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1','s2'},dom,[1 0],struct('in',1,'out',1),struct('psatz',1));
check_Zop(Zop,nrows_3pi(1,false(1,3)),{{(0:1)'}},{'s1'});
g = {@(s) (b-s(1)), @(s) (s(1)-a)};
spec = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,1-1)},'g',{g});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(~ismember('s2',Pmat{1}.varname),'Pmat must not depend on the unused variable s2.');
end

function case_unsorted_variables()
% Variables listed in non-alphabetical order, with their own domains: the
% psatz multipliers and monomials of Pmat must follow the variable names.
% The dummy degree is 0: for a positive dummy degree, monopmat (sopvar) does
% not return the same operator as for the alphabetical order (it warns
% 'sopvar:noncanonicalMultiplier'); Zop is not checked here.
rng(13);
dom = [-1 0.7; 0.5 2];          % s2 in [-1,0.7], s1 in [0.5,2]
vn = {'s2','s1'};
prog0 = new_prog(vn);
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,vn,dom,[1 1],struct('in',[0 0],'out',[2 1]),struct('psatz',[1 0]));
g = {@(s) (dom(1,2)-s(1)), @(s) (s(1)-dom(1,1))};      % s(1) is the value of s2
spec = struct('varnames',{vn},'dom',dom,'E',{uniform_E(1,2,[2-1,1])},'g',{g});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_positive_form()
% <Zop*x,Pmat*Zop*x>_{L2}, computed by polynomial arithmetic for arbitrary
% polynomial x (apply_op, then integration), must equal the sum of
% squares
%   sum_k int g_k(s)*|Rk*(I x Z(s))*(Zop*x)(s)|^2 ds,     Qk = Rk'*Rk,
% obtained from the Cholesky factors Rk of the Gram matrices and the
% monomials Z written out by hand, without using Pmat. In particular, it is
% positive: Pop = Zop'*Pmat*Zop is a positive operator.
rng(14);
a = 0.5;    b = 2;
% 1D with 2 states (second vector-valued), 1D with psatz, 2D, excluded term
cs = {};
cs{end+1} = struct('dims',[1;2],'vn',{{'s1'}},'dom',[a b],'dep',[1;1],...
    'degs',struct('in',{{1,2}},'out',{{2,1}}),'opts',struct(),'g',{{@(s) 1}});
cs{end+1} = struct('dims',[1;2],'vn',{{'s1'}},'dom',[a b],'dep',[1;1],...
    'degs',struct('in',{{1,1}},'out',{{3,2}}),'opts',struct('psatz',1),...
    'g',{{@(s) b-s(1), @(s) s(1)-a}});
cs{end+1} = struct('dims',1,'vn',{{'s1','s2'}},'dom',[a b;-1 0.7],'dep',[1 1],...
    'degs',struct('in',[1 1],'out',[1 0]),'opts',struct(),'g',{{@(s) 1}});
cs{end+1} = struct('dims',1,'vn',{{'s1'}},'dom',[a b],'dep',1,...
    'degs',struct('in',1,'out',1),'opts',struct('exclude',[true;false;false]),'g',{{@(s) 1}});
% Mixed dependence: states in their own variables, Pmat in all used ones.
% ip: also compute the value through innerprod(Zop*Zx,Pmat), the route of the
% programs (scalar states, <=2 variables, no excluded terms)
dm2 = [a b;-1 0.7];     dm3 = [a b;-1 0.7;1 3];
gn = {@(s) 1};
g01 = {@(s) dm2(2,2)-s(2), @(s) s(2)-dm2(2,1)};                            % psatz [0 1] in s2: (b2-s2), (s2-a2)
cs{end+1} = struct('dims',[1;1],'vn',{{'s1','s2'}},'dom',dm2,'dep',[1 0;1 1],...
    'degs',struct('in',[1 1],'out',[1 1]),'opts',struct(),'g',{gn},'ip',true);
cs{end+1} = struct('dims',[1;1],'vn',{{'s1','s2'}},'dom',dm2,'dep',[1 0;1 1],...
    'degs',struct('in',[1 1],'out',[1 2]),'opts',struct('psatz',[0 1]),'g',{g01},'ip',true);
cs{end+1} = struct('dims',[1;1],'vn',{{'s1','s2'}},'dom',dm2,'dep',[1 0;0 1],...
    'degs',struct('in',[1 1],'out',[1 1]),'opts',struct(),'g',{gn},'ip',true);
cs{end+1} = struct('dims',[1;1;1],'vn',{{'s1','s2'}},'dom',dm2,'dep',[1 0;0 1;1 1],...
    'degs',struct('in',[1 0],'out',[1 2]),'opts',struct('psatz',[0 1]),'g',{g01},'ip',true);
cs{end+1} = struct('dims',[2;3],'vn',{{'s1','s2'}},'dom',dm2,'dep',[1 0;1 1],...
    'degs',struct('in',[1 0],'out',[1 1]),'opts',struct(),'g',{gn});
cs{end+1} = struct('dims',[1;2],'vn',{{'s1','s2'}},'dom',dm2,'dep',[0 0;1 1],...
    'degs',struct('in',[1 1],'out',[1 0]),'opts',struct(),'g',{gn});
cs{end+1} = struct('dims',[1;1],'vn',{{'s1','s2','s3'}},'dom',dm3,'dep',[1 0 0;1 1 0],...
    'degs',struct('in',[1 1 1],'out',[1 1 1]),'opts',struct(),'g',{gn});
cs{end+1} = struct('dims',[1;1],'vn',{{'s1','s2','s3'}},'dom',dm3,'dep',[1 1 0;1 1 1],...
    'degs',struct('in',[1 0 1],'out',[0 0 0]),'opts',struct(),'g',{gn});
cs{end+1} = struct('dims',[1;1;1],'vn',{{'s1','s2','s3'}},'dom',dm3,'dep',[1 0 0;0 1 1;1 0 1],...
    'degs',struct('in',[1 0 1],'out',[0 0 0]),'opts',struct(),'g',{gn});
for c=1:numel(cs)
    cc = cs{c};
    prog0 = new_prog(cc.vn);
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,cc.dims,cc.vn,cc.dom,cc.dep,cc.degs,cc.opts);
    q = numel(cc.dims);     N = numel(cc.vn);       T = numel(cc.g);
    used = any(cc.dep>0,1);
    m = Zop.dims.out(:);
    assert(prog.var.num==T,'Case %d: expected %d declared SOS variables.',c,T);
    % Monomials Z: those of monopmat's interpretation of degs.out, cut in
    % degree by ceil(order/2) if there is a psatz multiplier
    E = ref_exponents(cc.dims,cc.vn,cc.dom,cc.dep,cc.degs);
    if isfield(cc.opts,'psatz')
        ps = cc.opts.psatz.*ones(1,N);      % order in each variable
        ps = ps(used);                      % ... on which a state depends
        for i=1:q
            for k=1:nnz(used)
                E{i}{k} = E{i}{k}(E{i}{k}<=max(E{i}{k})-ceil(ps(k)/2));
            end
        end
    end
    Zs = cell(q,1);
    for i=1:q
        Zs{i} = zvec_poly(E{i},cc.vn(used));
    end
    % Names of the decision variables in each Gram matrix, from direct
    % sosquadvar calls in the order of declaration
    pr = prog0;     Qn = cell(1,T);
    for t=1:T
        [pr,~,Qn{t}] = sosquadvar(pr,Zs,Zs,m,m,'pos');
    end
    nz = cellfun(@(e) prod(cellfun(@numel,e)),E);
    off = [0;cumsum(m(:).*nz(:))];
    % Random PSD Gram matrices Qk = Rk'*Rk, assigned to the decision variables
    names = cell(0,1);      vals = zeros(0,1);
    R = cell(1,T);
    for t=1:T
        Qt = random_psd(off(end));
        R{t} = chol(Qt);
        for i=1:q
            for j=1:q
                nm = Qn{t}{i,j};
                vl = Qt(off(i)+1:off(i+1),off(j)+1:off(j+1));
                names = [names;nm(:)];  %#ok<AGROW>
                vals = [vals;vl(:)];    %#ok<AGROW>
            end
        end
    end
    [names,iu] = unique(names,'stable');
    d0map = containers.Map(names,num2cell(vals(iu)));
    % Pmat at the assigned Gram matrices (once: eval_dpvar_at is slow), and
    % the functional innerprod(Zop*Zx,Pmat) (once: independent of the states)
    Pe = cell(q,q);
    for i=1:q
        for j=1:q
            Pe{i,j} = eval_dpvar_at(Pmat{i,j},d0map);
        end
    end
    rows = cell(q,1);
    for i=1:q
        rows{i} = horzcat(Pe{i,:});
    end
    Pfull = vertcat(rows{:});
    ipf = [];
    if isfield(cc,'ip') && cc.ip
        ipf = innerprod_functional(Pmat,Zop,cc);
        ipf.d0map = d0map;
    end
    for r=1:3
        [V,Vexp,Vip] = quad_form(Pfull,Zop,ipf,cc,r,Zs,R,m);
        assert(V>0,'Case %d: <Zop x,Pmat Zop x> must be positive; got %g.',c,V);
        assert(abs(V-Vexp)<=1e-8*(1+abs(Vexp)),...
            'Case %d: <Zop x,Pmat Zop x>=%g differs from the sum of squares %g.',c,V,Vexp);
        if isfield(cc,'ip') && cc.ip
            assert(abs(Vip-V)<=1e-8*(1+abs(V)),...
                'Case %d: innerprod(Zop*Zx,Pmat)=%g differs from the polynomial value %g.',c,Vip,V);
        end
    end
end
end

function [V,Vexp,Vip] = quad_form(Pfull,Zop,ipf,cc,seed,Zs,R,m)
% V = int <Zop*x,Pmat*Zop*x> for a random polynomial x, by polynomial
% arithmetic, with Pfull = Pmat at the assigned Gram matrices, and the sum of
% squares Vexp it must equal. State x_i is a polynomial in its own variables;
% the integral is over the variables on which some state depends. Vip: V
% through innerprod(Zop*Zx,Pmat) (NaN if ipf is empty).
rs = RandStream('mt19937ar','Seed',100+seed);
q = numel(cc.dims);     N = numel(cc.vn);
used = find(any(cc.dep>0,1));
svars = polynomial(cc.vn(:));
X = cell(1,q);
for i=1:q
    own = find(cc.dep(i,:)>0);
    comps = cell(cc.dims(i),1);
    for r=1:cc.dims(i)
        p = polynomial(randn(rs)+2);
        for k=own
            p = p + randn(rs)*svars(k) + randn(rs)*svars(k)^2;
        end
        if numel(own)>1, p = p + randn(rs)*svars(own(1))*svars(own(2)); end
        comps{r} = p;
    end
    X{i} = vertcat(comps{:});
end
Y = apply_op(fix_empty_params_zop(Zop),X);
% Using Pmat (Pfull: the blocks at the assigned Gram matrices)
F = Y.'*Pfull*Y;
for k=used
    F = int(F,svars(k),cc.dom(k,1),cc.dom(k,2));
end
V = double(F);
Vip = NaN;
if ~isempty(ipf)
    Vip = innerprod_eval(ipf,X);
end
% Using the Cholesky factors: (I x Z)*y = kron(y,Z) for each block
offm = [0;cumsum(m(:))];
blocks = cell(q,1);
for i=1:q
    yi = Y(offm(i)+1:offm(i+1));
    parts = cell(m(i),1);
    for r=1:m(i)
        parts{r} = yi(r)*Zs{i};
    end
    blocks{i} = vertcat(parts{:});
end
Wv = vertcat(blocks{:});
Vexp = 0;
for t=1:numel(cc.g)
    Wt = R{t}*Wv;
    Ft = cc.g{t}(svars)*(Wt.'*Wt);
    for k=used
        Ft = int(Ft,svars(k),cc.dom(k,1),cc.dom(k,2));
    end
    Vexp = Vexp + double(Ft);
end
end

function ipf = innerprod_functional(Pmat,Zop,cc)
% The functional <Zop*x,Pmat*Zop*x> through the route of the programs: the
% states as variables of a polyopvar Zx, u = Zop*Zx, and innerprod(u,Pmat,[]).
% Used for scalar states depending on <=2 variables: innerprod fails for
% mixed states in 3 variables (intop/mtimes_tensopvar).
q = numel(cc.dims);
Zx = [];
for i=1:q
    idx = find(cc.dep(i,:)>0);
    xi = polyopvar(sprintf('x%d',i),polynomial(cc.vn(idx)'),cc.dom(idx,:));
    if isempty(Zx), Zx = xi; else, Zx = [Zx;xi]; end    %#ok<AGROW>
end
u = Zop*Zx;
ipf = struct('V',innerprod(u,Pmat,[]),'dep_in',Zx.C.depmat.in,'vnames',{Zx.C.vars.in});
end

function Vip = innerprod_eval(ipf,X)
% Evaluate the functional of innerprod_functional at the states X (cell of
% polynomials) and the Gram matrices ipf.d0map
vals = struct();
for i=1:numel(X)
    vals.(sprintf('x%d',i)) = X{i};
end
Vd = eval_V(ipf.V,vals,ipf.dep_in,ipf.vnames);
Vip = double(eval_dpvar_at(Vd,ipf.d0map));
end

function Y = eval_V(V,vals,dep_in,vnames)
% Evaluate the quadratic functional V(x) (polyopvar, from innerprod) at the
% states vals.(name), polynomials in the variables vnames of the state. The
% dummy variable of variable s_k in the jth state factor depending on s_k is
% named s_k_j (as in apply_op), and by s_k itself for a single factor.
Y = 0;
nZ = size(V.Z,1);
Zf = full(V.Z);
for k=1:nZ
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    if d==0
        Y = Y + apply_op(op,{1});
    elseif isa(op,'intop')
        slot = zeros(1,numel(vnames));
        v = 1;
        for j=1:d
            xj = vals.(V.varname{idcs(j)});
            for mm=1:numel(vnames)
                if dep_in(idcs(j),mm)
                    slot(mm) = slot(mm)+1;
                    if d>1
                        isv = strcmp(xj.varname,vnames{mm});
                        xj.varname(isv) = {[vnames{mm},'_',num2str(slot(mm))]};
                    end
                end
            end
            v = kron(v,xj);
        end
        Y = Y + apply_op(op,v);
    else
        xv = cell(1,d);
        for j=1:d, xv{j} = vals.(V.varname{idcs(j)}); end
        Y = Y + apply_op(op,xv);
    end
end
end

function case_sequential()
% Two declarations in the same program, as for a Lyapunov functional and
% its multiplier: the decision variables must not overlap, and each of the
% matrices is an SOS matrix in its own right
rng(16);
a = 0.5;    b = 2;
prog0 = new_prog({'s1'});
[prog1,P1,Z1] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',1,'out',2));
opts = struct('psatz',1,'exclude',[true;false;false]);
[prog2,P2,Z2] = piesos_poslpivar(prog1,1,{'s1'},[a b],1,struct('in',2,'out',3),opts);
assert(prog2.var.num==3,'Expected 1+2 declared SOS variables.');
assert(size(Z2,1)==nrows_3pi(2,[true false false]),'Wrong row count of the second Zop.');
spec1 = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,2)},'g',{no_multiplier()});
check_sos(prog1,P1,Z1,spec1,prog0);
spec2 = struct('varnames',{{'s1'}},'dom',[a b],'E',{uniform_E(1,1,3-1)},'g',{{@(s) b-s(1),@(s) s(1)-a}});
check_sos(prog2,P2,Z2,spec2,prog1);
assert(isempty(intersect(P1{1}.dvarname,P2{1}.dvarname)),'The two matrices must not share decision variables.');
end

function case_no_variables()
% No spatial variables: states in R^2 and R, and Pmat a constant PSD matrix
rng(17);
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[2;1],{},zeros(0,2),zeros(2,0));
assert(isequal(size(Pmat),[2,2]),'Pmat must be a 2 x 2 cell.');
assert(isequal(Zop.dims.out(:),[2;1]),'Zop must be the identity on the states.');
spec = struct('varnames',{cell(1,0)},'dom',zeros(0,2),'E',{{cell(1,0),cell(1,0)}},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_sdp_feasibility()
% Pmat as an SOS matrix variable in an actual SDP (sossolve, SeDuMi): the
% equality Pmat = Ptarget is feasible exactly if Ptarget can be written as
%   sum_k g_k*(I x Z)'*Qk*(I x Z),   Qk >= 0,
% with the psatz multipliers g_k, nonnegative on the domain only:
%  - a target built from a PSD Gram matrix is feasible, -I is not;
%  - g*I, g=(s-a)(b-s), is positive on [a,b] only: feasible with a psatz
%    multiplier, infeasible without, and -g*I is infeasible in both cases;
%  - for two states, the blocks of Pmat form ONE PSD matrix: the target
%    [I 2I; 2I I] has PSD diagonal blocks, but is infeasible;
%  - the same for products of multipliers in two variables.
rng(15);
a = 0.5;    b = 2;
s = polynomial({'s1'});
g = (s-a)*(b-s);
Z = [polynomial(1);s];
% (a) plain SOS in one variable
prog0 = new_prog({'s1'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',0,'out',1));
m = size(Zop,1);
IZ = kron_IZ(m,Z);
Q0 = random_psd(2*m);
assert(is_feasible(prog,{Pmat{1,1}-IZ.'*Q0*IZ}),'A target built from a Gram matrix must be feasible.');
assert(~is_feasible(prog,{Pmat{1,1}+eye(m)}),'The target -I is not an SOS matrix, and must be infeasible.');
% (b) psatz multipliers
[progP,PmatP,ZopP] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',0,'out',2),struct('psatz',2));
mP = size(ZopP,1);
assert(is_feasible(progP,{PmatP{1,1}-g*eye(mP)}),'g*I must be feasible with a psatz multiplier.');
assert(~is_feasible(progP,{PmatP{1,1}+g*eye(mP)}),'-g*I must be infeasible with a psatz multiplier.');
[progN,PmatN,ZopN] = piesos_poslpivar(prog0,1,{'s1'},[a b],1,struct('in',0,'out',2));
assert(~is_feasible(progN,{PmatN{1,1}-g*eye(size(ZopN,1))}),'g*I is not positive on the real line: infeasible without psatz.');
% (c) two states, one joint PSD constraint
[prog2,Pmat2,Zop2] = piesos_poslpivar(prog0,[1;1],{'s1'},[a b],[1;1],struct('in',0,'out',1));
m2 = Zop2.dims.out(:);
IZ2 = blkdiag_poly(kron_IZ(m2(1),Z),kron_IZ(m2(2),Z));
Pt = IZ2.'*random_psd(2*sum(m2))*IZ2;
r1 = 1:m2(1);       r2 = m2(1)+(1:m2(2));
assert(is_feasible(prog2,{Pmat2{1,1}-Pt(r1,r1),Pmat2{1,2}-Pt(r1,r2),Pmat2{2,1}-Pt(r2,r1),Pmat2{2,2}-Pt(r2,r2)}),...
    'A joint Gram matrix target for two states must be feasible.');
T = kron([1 2;2 1],eye(m2(1)));
assert(~is_feasible(prog2,{Pmat2{1,1}-T(r1,r1),Pmat2{1,2}-T(r1,r2),Pmat2{2,1}-T(r2,r1),Pmat2{2,2}-T(r2,r2)}),...
    'The blocks of Pmat must form one PSD matrix: [I 2I;2I I] is infeasible.');
% (d) two variables, psatz [2 2]: g1*g2*I positive on the rectangle only
dom = [a b;-1 0.7];
x2 = polyopvar('x1',polynomial({'s1','s2'}),dom);
progD = piesos_program(x2);
[progD,PmatD,ZopD] = piesos_poslpivar(progD,1,{'s1','s2'},dom,[1 1],struct('in',[0 0],'out',[2 2]),struct('psatz',[2 2]));
mD = size(ZopD,1);
sv = polynomial({'s1','s2'});
g2 = (sv(1)-dom(1,1))*(dom(1,2)-sv(1))*(sv(2)-dom(2,1))*(dom(2,2)-sv(2));
assert(is_feasible(progD,{PmatD{1,1}-g2*eye(mD)}),'g1*g2*I must be feasible with psatz [2 2].');
assert(~is_feasible(progD,{PmatD{1,1}+g2*eye(mD)}),'-g1*g2*I must be infeasible with psatz [2 2].');
end

function case_errors()
prog = new_prog({'s1'});
assert_err(@() piesos_poslpivar(prog,1,struct('in',{{'s1'}},'out',{{'s2'}}),[0 1],1,1,struct()),...
    'must match');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},struct('in',[0 1],'out',[0 2]),1,1,struct()),...
    'must match');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],struct('in',1,'out',2),1,struct()),...
    'must match');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('psatz',1,'otype',[2 2])),...
    'currently not supported');
% Wrong sizes are rejected by monopmat
assert_err(@() piesos_poslpivar(prog,[1;1],{'s1'},[0 1],1,1,struct()),'');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1;0 1],1,1,struct()),'');
% Invalid degrees
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,struct('in',1,'out',-1),struct()),'nonnegative integers');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,struct('in',1,'out',1.5),struct()),'nonnegative integers');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,struct('in',1),struct()),'fields ''in'' and ''out''');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,struct('in',1,'out',[1 2]),struct()),'1 x M');
% Invalid psatz
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('psatz',-1)),'nonnegative integers');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('psatz',1.5)),'nonnegative integers');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('psatz',[1 2])),'1 x N array');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('psatz','a')),'1 x N array');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('psatz',[false;false])),'must be true');
prog2 = new_prog({'s1','s2'});
assert_err(@() piesos_poslpivar(prog2,1,{'s1','s2'},[0 1;0 1],[1 0],1,struct('psatz',[1 1])),'no state depends');
assert_err(@() piesos_poslpivar(prog2,1,{'s1','s2'},[0 1;0 1],[1 1],1,struct('psatz',true(2,2,2))),'1 x N array');
% States depending on different variables are accepted: see case_mixed_*
piesos_poslpivar(prog2,[1;1],{'s1','s2'},[0 1;0 1],[1 0;1 1],1,struct());
piesos_poslpivar(prog2,[1;1],{'s1','s2'},[0 1;0 1],[1 1;0 0],1,struct());
% Operators on the left and right with equal dimensions, that differ: exclude the
% lower vs the upper integral, or dummy monomials 1,s vs 1,s^2 (same row count)
exL = struct('exclude',[false;true;false]);
exR = struct('exclude',[false;false;true]);
dL = struct('in',{{{[0;1]}}},'out',{{{[0;1]}}});
dR = struct('in',{{{[0;2]}}},'out',{{{[0;1]}}});
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,1),'must be specified as a');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('L',exL,'R',exR)),'same operators');
assert_err(@() piesos_poslpivar(prog,1,{'s1'},[0 1],1,struct('L',dL,'R',dR),struct()),'same operators');
% Equal operators are accepted, and give the same result as common options
[progA,PmatA,ZopA] = piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,struct('L',exL,'R',exL));
[progB,PmatB,ZopB] = piesos_poslpivar(prog,1,{'s1'},[0 1],1,1,exL);
assert(isequal(progA,progB) && isequal(ZopA,ZopB) && isequal(PmatA{1}.C,PmatB{1}.C),...
    'opts.L=opts.R must give the same result as common options.');
[progA,PmatA,ZopA] = piesos_poslpivar(prog,1,{'s1'},[0 1],1,struct('L',dL,'R',dL),struct());
[progB,PmatB,ZopB] = piesos_poslpivar(prog,1,{'s1'},[0 1],1,dL,struct());
assert(isequal(progA,progB) && isequal(ZopA,ZopB) && isequal(PmatA{1}.C,PmatB{1}.C),...
    'degs.L=degs.R must give the same result as common degrees.');
end


function [vn,dom] = mixed_vars(N)
% Variables s1..sN with their own domains
vn = arrayfun(@(k) sprintf('s%d',k),1:N,'UniformOutput',false);
dom = [0.5 2; -1 0.7; 1 3];
dom = dom(1:N,:);
end

function m = nrows_nd(r,excl,N)
% Rows of the basis of a state depending on variables with dummy degrees r (one
% per variable, in a list of N variables): the product over the variables of
% the state, less the excluded terms. A 3-vector excl = logical
% [mult,lower,upper] excludes the terms with an excluded kind in ANY variable
% of the state (OR), as monopmat applies it (10/21/2026; was AND over the list,
% with the term a multiplier in the variables the state does not depend on):
% the rows are the product over the variables of the state of the rows of the
% kinds not excluded. A state without variables has its multiplier only.
if isempty(r)
    m = ~excl(1)*1;
else
    m = prod(arrayfun(@(rk) nrows_3pi(rk,excl),r));
end
end

function case_mixed_structure()
% States depending on different variables: Zop has, per state, the rows of its
% own input variables (product of nrows_3pi), ZL=0, dummy monomials (0:r)' in
% the state's variables only, and output variables covering ALL used variables
% (the common space). Every Pmat{i,j} is a function of all used variables with
% monomials (0:d)'.
rng(18);
cs = {};
%                dims        depmat                            dummy r    primary d
cs{end+1} = {[1;1],      [1 0;1 1],                        [1 1 0],   [1 1 0]};     % nested
cs{end+1} = {[1;2],      [1 0;0 1],                        [1 2 0],   [2 1 0]};     % disjoint
cs{end+1} = {[1;1],      [0 0;1 1],                        [1 1 0],   [1 1 0]};     % ghost state
cs{end+1} = {[1;2;1],    [0 0;1 0;1 1],                    [1 1 0],   [1 0 0]};     % ghost, 1D, 2D
cs{end+1} = {[2;3],      [1 0;1 1],                        [1 0 0],   [1 0 0]};     % vector-valued
cs{end+1} = {[1;1],      [1 1 0;1 1 1],                    [1 0 1],   [1 0 0]};     % 3 variables
cs{end+1} = {[1;1;1],    [1 0 0;0 1 1;1 0 1],              [1 0 1],   [1 0 0]};     % 3 variables
cs{end+1} = {[1;1],      [1 0 0;1 1 0],                    [1 1 1],   [1 1 1]};     % s3 unused
cs{end+1} = {[1;1],      [0 1;1 0],                        [1 1 0],   [0 2 0]};     % swapped
for c=1:numel(cs)
    [dims,dep,r,d] = cs{c}{:};
    q = numel(dims);    N = size(dep,2);
    [vn,dom] = mixed_vars(N);
    used = any(dep>0,1);
    prog0 = new_prog(vn);
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,dims,vn,dom,dep,struct('in',r(1:N),'out',d(1:N)),struct());
    m = zeros(q,1);     ZR = cell(q,1);
    for i=1:q
        own = find(dep(i,:)>0);
        m(i) = dims(i)*nrows_nd(r(own),false(1,3),N);
        ZR{i} = cell(1,N);
        for k=own, ZR{i}{k} = (0:r(k))'; end
    end
    check_Zop(Zop,m,ZR,vn,vn(used));
    spec = struct('varnames',{vn(used)},'dom',dom(used,:),'E',{uniform_E(q,nnz(used),d(used))},'g',{no_multiplier()});
    check_sos(prog,Pmat,Zop,spec,prog0);
    % Every block of Pmat depends on all used variables where d>0 (also if the
    % state does not), and on no unused variable
    for i=1:q
        for j=1:q
            assert(all(ismember(Pmat{i,j}.varname,vn(used))),'Case %d: Pmat{%d,%d} depends on an unused variable.',c,i,j);
        end
    end
    for k=find(used & d(1:N)>0)
        for i=1:q
            assert(any(strcmp(Pmat{i,i}.varname,vn{k})),...
                'Case %d: Pmat{%d,%d} must depend on s%d (common space).',c,i,i,k);
        end
    end
end
% depmat=2 for one state (two factors), 1 for the other: the 5x5 first block
% and the product of the factors' rows for the second
dom = [0.5 2;-1 0.7];
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],{'s1','s2'},dom,[1 0;2 2],struct('in',[1 1],'out',[1 0]),struct());
m2 = (nrows_3pi(1,false(1,3))^2)*(nrows_3pi(0,false(1,3))^2);
assert(isequal(Zop.dims.out(:),[nrows_3pi(1,false(1,3));m2]),'depmat [1 0;2 2]: wrong row counts %s.',mat2str(Zop.dims.out(:)'));
assert(numel(Zop.ops{1,1}.ops)==1 && numel(Zop.ops{2,2}.ops)==2,'Expected 1 and 2 factors.');
for i=1:2
    for f=1:numel(Zop.ops{i,i}.ops)
        assert(all(cellfun(@(z) isequal(z(:),0),Zop.ops{i,i}.ops{f}.ZL)),'ZL must be 0 (block %d, factor %d).',i,f);
    end
end
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(2,2,[1 0])},'g',{no_multiplier()});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_mixed_exclude_sep()
% Excluded terms and separability with states of different dependence (see
% nrows_nd): e.g. for the state of x1(s1) in the list (s1,s2), exclude [0 1 0]
% removes the lower integral, as monopmat applies the 3-vector along the
% variables of each state (10/21/2026; was: left it in).
rng(19);
dom = [0.5 2;-1 0.7];
r = [2 1];
for dep = {[1 0;1 1],[0 0;1 1],[1 0;0 1]}
    dep = dep{1};                                                       %#ok<FXSET>
    q = size(dep,1);
    for excl = {[1 0 0],[0 1 0],[0 0 1],[0 1 1],[1 1 0]}
        excl = logical(excl{1});                                        %#ok<FXSET>
        prog0 = new_prog({'s1','s2'});
        [prog,Pmat,Zop] = piesos_poslpivar(prog0,ones(q,1),{'s1','s2'},dom,dep,struct('in',r,'out',[1 1]),...
            struct('exclude',excl(:)));
        m = zeros(q,1);
        for i=1:q
            m(i) = nrows_nd(r(dep(i,:)>0),excl,2);
        end
        assert(isequal(Zop.dims.out(:),m),'dep=%s, exclude=%s: rows %s expected, got %s.',mat2str(dep),...
            mat2str(excl),mat2str(m'),mat2str(Zop.dims.out(:)'));
        if any(m==0)
            continue        % the multiplier of a state without variables is its only term
        end
        spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(q,2,[1 1])},'g',{no_multiplier()});
        check_sos(prog,Pmat,Zop,spec,prog0);
    end
    % separable: the upper integral shares the parameter of the lower one,
    % in each variable of the state, also for a state that depends on s2 only
    prog0 = new_prog({'s1','s2'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,ones(q,1),{'s1','s2'},dom,dep,struct('in',r,'out',[1 1]),struct('sep',true));
    m = zeros(q,1);
    for i=1:q
        m(i) = prod(arrayfun(@(rk) nrows_3pi(rk,[false false true]),r(dep(i,:)>0)));
    end
    assert(isequal(Zop.dims.out(:),m),'dep=%s, sep: rows %s expected, got %s.',mat2str(dep),mat2str(m'),mat2str(Zop.dims.out(:)'));
    spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(q,2,[1 1])},'g',{no_multiplier()});
    check_sos(prog,Pmat,Zop,spec,prog0);
end
end

function g = psatz_terms(dom,ps,lin)
% Multipliers prod_k (s_k-a_k)^jk*(b_k-s_k)^(qk-jk) for the terms with linear
% index lin (first variable fastest) in the array of size ps+1; ps(k)=order
% in variable k, as functions of the vector s of the variables
g = cell(1,numel(lin));
for t=1:numel(lin)
    jj = cell(1,numel(ps));
    [jj{:}] = ind2sub(ps+1,lin(t));
    j = cell2mat(jj)-1;
    g{t} = @(s) prod((s(:)'-dom(:,1)').^j.*(dom(:,2)'-s(:)').^(ps-j));
end
end

function case_mixed_psatz()
% psatz multipliers in variables on which some of the states do not depend:
% the multiplier acts on every block of Pmat, as all blocks are functions of
% the common variables. Monomials cut by ceil(qk/2) in the variable k.
rng(20);
dom = [0.5 2;-1 0.7];
% (a) multiplier in s2 only, on the state that depends on s1 only and on both
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],{'s1','s2'},dom,[1 0;1 1],struct('in',[1 0],'out',[1 2]),struct('psatz',[0 1]));
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(2,2,[1,2-1])},...
    'g',{psatz_terms(dom,[0 1],1:2)});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(prog.var.num==2,'Expected 2 SOS variables.');
assert(any(strcmp(Pmat{1,1}.varname,'s2')),'Pmat{1,1} must depend on s2, on which state 1 does not depend.');
% (b) orders [1 2] in both variables, 6 terms
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],{'s1','s2'},dom,[1 0;1 1],struct('in',[1 0],'out',[2 3]),struct('psatz',[1 2]));
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(2,2,[2-1,3-1])},'g',{psatz_terms(dom,[1 2],1:6)});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(prog.var.num==6,'Expected 6 SOS variables.');
% (c) logical array, 2 x 3: selected terms only
L = false(2,3);     L(1,1) = true;      L(2,2) = true;      L(1,3) = true;
prog0 = new_prog({'s1','s2'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1;1],{'s1','s2'},dom,[1 0;0 1;1 1],struct('in',[0 0],'out',[2 2]),struct('psatz',L));
spec = struct('varnames',{{'s1','s2'}},'dom',dom,'E',{uniform_E(3,2,[2-1,2-1])},'g',{psatz_terms(dom,[1 2],find(L)')});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(prog.var.num==3,'Expected one SOS variable for each selected multiplier.');
% (d) scalar psatz with an unused variable: applied to the used variables only
dom3 = [0.5 2;-1 0.7;1 3];
prog0 = new_prog({'s1','s2','s3'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],{'s1','s2','s3'},dom3,[1 0 0;1 1 0],struct('in',[1 1 1],'out',[2 2 0]),struct('psatz',1));
spec = struct('varnames',{{'s1','s2'}},'dom',dom3(1:2,:),'E',{uniform_E(2,2,[2-1,2-1])},'g',{psatz_terms(dom3(1:2,:),[1 1],1:4)});
check_sos(prog,Pmat,Zop,spec,prog0);
assert(prog.var.num==4 && ~ismember('s3',Pmat{1}.varname),'Expected 4 SOS variables and no dependence on s3.');
% (e) 3 variables, multiplier in s3 on which the first state does not depend
prog0 = new_prog({'s1','s2','s3'});
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],{'s1','s2','s3'},dom3,[1 1 0;1 1 1],struct('in',[1 0 0],'out',[0 0 1]),struct('psatz',[0 0 1]));
spec = struct('varnames',{{'s1','s2','s3'}},'dom',dom3,'E',{uniform_E(2,3,[0 0 1-1])},'g',{psatz_terms(dom3,[0 0 1],1:2)});
check_sos(prog,Pmat,Zop,spec,prog0);
% (f) variables in non-alphabetical order, mixed (dummy degree 0, see case_unsorted_variables)
vn = {'s2','s1'};   dm = dom([2 1],:);
prog0 = new_prog(vn);
[prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],vn,dm,[0 1;1 1],struct('in',[0 0],'out',[2 1]),struct('psatz',[1 0]));
spec = struct('varnames',{vn},'dom',dm,'E',{uniform_E(2,2,[2-1,1])},'g',{psatz_terms(dm,[1 0],1:2)});
check_sos(prog,Pmat,Zop,spec,prog0);
end

function case_mixed_sdp()
% Mixed dependence in SDPs, Zop of dummy degree 0 for states x1(s1) and
% x2(s1,s2), common output space (s1,s2):
%  - a target from ONE Gram matrix over both blocks (of different size) is
%    feasible, a target PSD on each diagonal block but not jointly
%    ([I 2E;2E' I], E=[I 0]) is not;
%  - g(s2)=(s2-a2)(b2-s2), positive on [a2,b2] only, in a variable on which x1
%    does not depend: g(s2)*I on ALL blocks is feasible with psatz [0 2],
%    infeasible without, and -g(s2)*I is infeasible with it;
%  - the same with the multiplier in s1 (state 1 and state 2 both depend on
%    s1) as a control, and for a disjoint pair x1(s1), x2(s2).
rng(21);
dom = [0.5 2;-1 0.7];
sv = polynomial({'s1','s2'});
for dep = {[1 0;1 1],[1 0;0 1]}
    dep = dep{1};                                                       %#ok<FXSET>
    % (a) one Gram matrix over both blocks
    prog0 = new_prog({'s1','s2'});
    [prog,Pmat,Zop] = piesos_poslpivar(prog0,[1;1],{'s1','s2'},dom,dep,struct('in',[0 0],'out',[1 1]),struct());
    m = Zop.dims.out(:);
    E = {(0:1)',(0:1)'};
    Z = zvec_poly(E,{'s1','s2'});
    IZ = blkdiag_poly(kron_IZ(m(1),Z),kron_IZ(m(2),Z));
    Pt = IZ.'*random_psd(numel(Z)*sum(m))*IZ;
    r1 = 1:m(1);        r2 = m(1)+(1:m(2));
    assert(is_feasible(prog,{Pmat{1,1}-Pt(r1,r1),Pmat{1,2}-Pt(r1,r2),Pmat{2,1}-Pt(r2,r1),Pmat{2,2}-Pt(r2,r2)}),...
        'dep=%s: a target from one joint Gram matrix must be feasible.',mat2str(dep));
    Ep = [eye(m(1)),zeros(m(1),m(2)-m(1))];
    if m(2)<m(1), Ep = [eye(m(2));zeros(m(1)-m(2),m(2))]; end
    T = [eye(m(1)),2*Ep;2*Ep',eye(m(2))];
    assert(~is_feasible(prog,{Pmat{1,1}-T(r1,r1),Pmat{1,2}-T(r1,r2),Pmat{2,1}-T(r2,r1),Pmat{2,2}-T(r2,r2)}),...
        'dep=%s: blocks of Pmat must form one PSD matrix; [I 2E;2E'' I] is infeasible.',mat2str(dep));
    % (b) psatz multiplier in s2: target g(s2)*I on all blocks
    g2 = (sv(2)-dom(2,1))*(dom(2,2)-sv(2));
    g1 = (sv(1)-dom(1,1))*(dom(1,2)-sv(1));
    for kvar=[2 1]
        gk = g2;    ps = [0 2];
        if kvar==1, gk = g1;    ps = [2 0]; end
        prog0 = new_prog({'s1','s2'});
        [progP,PmatP,ZopP] = piesos_poslpivar(prog0,[1;1],{'s1','s2'},dom,dep,struct('in',[0 0],'out',[2 2]),struct('psatz',ps));
        mm = ZopP.dims.out(:);
        rr = {1:mm(1),mm(1)+(1:mm(2))};
        Gt = gk*eye(sum(mm));
        feas = @(sgn) is_feasible(progP,[cellfun(@(i,j) PmatP{i,j}-sgn*Gt(rr{i},rr{j}),...
            num2cell(repelem(1:2,2)),num2cell(repmat(1:2,1,2)),'UniformOutput',false)]);
        assert(feas(1),'dep=%s: g(s%d)*I must be feasible with psatz %s.',mat2str(dep),kvar,mat2str(ps));
        assert(~feas(-1),'dep=%s: -g(s%d)*I must be infeasible with psatz %s.',mat2str(dep),kvar,mat2str(ps));
        [progN,PmatN,ZopN] = piesos_poslpivar(prog0,[1;1],{'s1','s2'},dom,dep,struct('in',[0 0],'out',[2 2]),struct());
        mn = ZopN.dims.out(:);
        rn = {1:mn(1),mn(1)+(1:mn(2))};
        assert(~is_feasible(progN,cellfun(@(i,j) PmatN{i,j}-Gt(rn{i},rn{j}),...
            num2cell(repelem(1:2,2)),num2cell(repmat(1:2,1,2)),'UniformOutput',false)),...
            'dep=%s: g(s%d)*I is not positive on the real line, so infeasible without psatz.',mat2str(dep),kvar);
    end
end
end


%% ---- helpers for the SDP case ------------------------------------------
function feas = is_feasible(prog,exprs)
% Feasibility of the program with the equalities exprs{k}==0 added
for k=1:numel(exprs)
    prog = soseq(prog,exprs{k});
end
out = evalc('[prog,info] = sossolve(prog);');       %#ok<ASGLU>
feas = (info.pinf==0 && info.dinf==0);
end

function IZ = kron_IZ(m,Z)
% The (m*nz) x m polynomial matrix I_m x Z, with Z an nz x 1 polynomial
nz = size(Z,1);
cols = cell(1,m);
for r=1:m
    cols{r} = vertcat(polynomial(zeros(nz*(r-1),1)),Z,polynomial(zeros(nz*(m-r),1)));
end
IZ = horzcat(cols{:});
end

function C = blkdiag_poly(A,B)
% Block-diagonal concatenation of polynomial matrices
C = [A,polynomial(zeros(size(A,1),size(B,2)));polynomial(zeros(size(B,1),size(A,2))),B];
end

%% ---- workaround for apply_sopvar (as in Test_piesos_sosvar.m) ---------
function Zop2 = fix_empty_params_zop(Zop)
% monopmat leaves params{k}=[] for any excluded multiplier/lower/upper
% combo, which apply_sopvar cannot handle: expand to zero sparse matrices
% (only needed for the ground truth evaluation)
Zop2 = Zop;
[p,d] = size(Zop.ops);
for i=1:p
    for j=1:d
        if ~isempty(Zop.ops{i,j})
            Zop2.ops{i,j} = fix_empty_params_tensopvar(Zop.ops{i,j});
        end
    end
end
end

function T2 = fix_empty_params_tensopvar(T)
T2 = T;
[nt,nf] = size(T.ops);
for i=1:nt
    for j=1:nf
        entry = T.ops{i,j};
        if isa(entry,'sopvar')
            T2.ops{i,j} = fill_empty_params_sopvar(entry);
        elseif isa(entry,'tensopvar')
            T2.ops{i,j} = fix_empty_params_tensopvar(entry);
        end
    end
end
end

function Pop2 = fill_empty_params_sopvar(Pop)
params = Pop.params;
nL = prod(cellfun(@numel,Pop.ZL));
nR = prod(cellfun(@numel,Pop.ZR));
nrow = Pop.dims(1)*nL;      ncol = Pop.dims(2)*nR;
for k=1:numel(params)
    if isempty(params{k})
        params{k} = sparse(nrow,ncol);
    end
end
Pop2 = sopvar(params,Pop.vars,Pop.ZL,Pop.ZR,Pop.dom,Pop.dims);
end
