function Test_piesos_quadvar()
% TEST_PIESOS_QUADVAR tests, using assertions, piesos_quadvar,
%   [prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog,ZL,ZR,opts),
% which declares V = <ZopL*ZL, Pmat*ZopR*ZR>, with ZopL, ZopR bases of
% tensor-PI operators acting on the distributed monomials in ZL, ZR, and
% Pmat an nL x nR cell of matrices in the primary variables (general,
% symmetric or SOS).
%
% Ground truths are independent of the function under test:
%  - V: <ZopL*zL, Pmat*ZopR*zR>_{L2} by polynomial arithmetic (apply_op on
%    ZopL, ZopR for random polynomial states, then integration over the
%    variables on which some monomial depends), with random values for the
%    decision variables of Pmat (looked up by name); compared with the
%    returned V at the same states and values, through its kernels
%    (apply_op, as in Test_piesos_sosvar). V is linear in the decision
%    variables, so random values test the full structure of Pmat and the
%    kernels;
%  - equal bases with opts.pos=true: program, Pmat, ZopL=ZopR=Zop and V of
%    piesos_sosvar (verified by Test_piesos_sosvar);
%  - the wiring to piesos_quadopvar: Pmat, ZopL, ZopR and the program must
%    be those of a direct call with the function spaces of ZL.C, ZR.C
%    written out by hand;
%  - the forms of opts (degrees in opts.degs or opts.L.degs, opts.R.degs,
%    other operator options) and the errors;
%  - the coefficients ZL.C, ZR.C of the bases (operators T, scalings, sums,
%    tensor products, vector selections): V is the quadratic form of the
%    entries of ZL, ZR, with Pmat indexed by their rows; the ground truth
%    builds the entries from the outputs of the coefficients; and the
%    unsupported coefficients are rejected.
% No programs are solved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_quadvar
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
% DJ (CC), 10/21/2026: Add case_basis_coefficients;
% DJ (CC), 10/22/2026: Replace case_basis_coefficients by
%                        case_coefficients, case_coefficient_errors: V for
%                        coefficients ZL.C, ZR.C, rows of Pmat (apply_side,
%                        eval_rows);
% DJ (CC), 10/22/2026: Add case_ode_state;
% DJ (CC), 10/23/2026: Bases from otimes (tensor products), not x1.*x1, which
%                        is rejected (new case in case_coefficient_errors);
% DJ (CC), 10/23/2026: dmonomials has the variables y, not the scalar rows:
%                        basis of the entries of y built by selection, new
%                        cases dmonomials(y,[1;2]), dmonomials(Tm*y,[1;2]);

% NOTES
% - piesos_sosvar wraps piesos_quadvar, so the comparison with it only
%   checks the forwarding of the inputs; V is verified by the direct
%   evaluation;
% - Ground truth needs fix_empty_params_zop (monopmat leaves params{k}=[]
%   for excluded multiplier/lower/upper combinations, which apply_op
%   cannot handle);
% - Random decision-variable values come from rng(seed) per check, so
%   reruns are reproducible.

tests ={@case_equals_sosvar, @case_general_equal_bases, @case_different_bases, ...
         @case_operators_differ, @case_vector_valued, @case_basis_forms, ...
         @case_degs_forms, @case_operator_options, @case_wiring, ...
         @case_coefficients, @case_ode_state, @case_coefficient_errors, @case_errors};
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
    error('Test_piesos_quadvar: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_quadvar: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function S = make_states()
% States in s1 in [0.5,2], s2 in [-1,0.7], and polynomial values at which to evaluate V
dom = [0.5 2;-1 0.7];
S.dom = dom;
S.x1 = polyopvar('x1',polynomial({'s1'}),dom(1,:));                 % x1(s1)
S.x2 = polyopvar('x2',polynomial({'s1','s2'}),dom);                 % x2(s1,s2)
S.x3 = polyopvar('x3',polynomial({'s2'}),dom(2,:));                 % x3(s2)
S.x4 = polyopvar('x4',polynomial({'s1'}),dom(1,:));                 % x4(s1)
S.y = polyopvar('y',2,polynomial({'s1'}),dom(1,:));                 % y(s1) in R^2
S.x5 = polyopvar('x5',polynomial({'s1','s2'}),dom);                 % x5(s1,s2)
s1 = polynomial({'s1'});    s1 = s1(1);
s2 = polynomial({'s2'});    s2 = s2(1);
S.X = struct('x1',1.5+0.7*s1-0.3*s1^2,...
             'x2',1.2+0.4*s1+0.9*s2-0.2*s1*s2+0.1*s2^2,...
             'x3',0.8-0.6*s2+0.5*s2^2,...
             'x4',-0.4+1.1*s1+0.2*s1^2,...
             'x5',0.9-0.3*s1+0.5*s2+0.4*s1*s2-0.1*s2^2,...
             'y',[1.1+0.3*s1-0.2*s1^2; 0.7-0.5*s1+0.4*s1^2]);
S.all = [S.x1;S.x2;S.x3;S.x4;S.x5;S.y];
end

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s" (expected "%s").',err.message,needle);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function m = Zdims(Zop)
m = Zop.dims.out(:);
end

function assert_same_output(A,B,label)
% Outputs A, B = {prog,V,Pmat,ZopL,ZopR} of two calls must be identical
assert(isequal(A{1},B{1}),'%s: programs differ.',label);
assert(isequal(A{2},B{2}),'%s: V differs.',label);
assert(isequal(size(A{3}),size(B{3})),'%s: sizes of Pmat differ.',label);
for k=1:numel(A{3})
    assert(isequal(A{3}{k}.C,B{3}{k}.C) && isequal(A{3}{k}.degmat,B{3}{k}.degmat) ...
        && isequal(A{3}{k}.varname,B{3}{k}.varname) && isequal(A{3}{k}.dvarname,B{3}{k}.dvarname),...
        '%s: Pmat{%d} differs.',label,k);
end
assert(isequal(A{4},B{4}),'%s: ZopL differs.',label);
assert(isequal(A{5},B{5}),'%s: ZopR differs.',label);
end

function d0map = random_decvars(Pmat)
% Random values for the decision variables in Pmat, by name
names = cell(0,1);
for k=1:numel(Pmat)
    names = [names;Pmat{k}.dvarname(:)];        %#ok<AGROW>
end
names = unique(names,'stable');
d0map = containers.Map(names,num2cell(0.5+rand(numel(names),1)));
end

function [Vtrue,Vlib] = check_V(V,Pmat,ZopL,ZopR,ZL,ZR,X,d0map,rowsL,rowsR)
% Direct evaluation of <ZopL*zL,Pmat*ZopR*zR> (polynomial arithmetic, no
% int/intop) and of V through its kernels, at the states X (struct: state
% name -> polynomial (vector)) and decision variable values d0map. The
% entries of zL, zR are those of the monomial bases ZL, ZR (identity
% coefficients), or, if rowsL, rowsR are given, the factors listed there, see
% eval_rows
if nargin<9 || isempty(rowsL),  YL = apply_side(ZopL,ZL,X);     else,   YL = eval_rows(ZopL,rowsL);     end
if nargin<10 || isempty(rowsR), YR = apply_side(ZopR,ZR,X);     else,   YR = eval_rows(ZopR,rowsR);     end
[q1,q2] = size(Pmat);
rows = cell(q1,1);
for i=1:q1
    blks = cell(1,q2);
    for j=1:q2
        blks{j} = eval_dpvar_at(Pmat{i,j},d0map);
    end
    rows{i} = horzcat(blks{:});
end
F = YL.'*vertcat(rows{:})*YR;
% integrate over the variables on which some monomial depends
[vn,dm] = used_vars(ZL,ZR);
for k=1:numel(vn)
    sk = polynomial(vn(k));
    F = int(F,sk(1),dm(k,1),dm(k,2));
end
Vtrue = double(F);
Vd = eval_V_generic(V,X,double(V.depmat.in)>0,V.vars.in);
Vlib = double(eval_dpvar_at(Vd,d0map));
end

function run_check(label,S,ZL,ZR,opts,seed,vs_sosvar,rowsL,rowsR)
% Declare V for the bases ZL, ZR (ZR=[]: same basis) and compare its value with
% the direct evaluation at random decision variable values. If vs_sosvar
% (ZL=ZR, opts.pos), also compare with piesos_sosvar (opts without pos). For
% coefficients other than the identity, the factors in each row of ZL, ZR are
% given in rowsL, rowsR, see eval_rows
if nargin<7,    vs_sosvar = false;      end
if nargin<8,    rowsL = [];     end
if nargin<9,    rowsR = [];     end
try
    run_check_body(label,S,ZL,ZR,opts,seed,vs_sosvar,rowsL,rowsR);
catch err
    fprintf('        FAILED: %s\n',label);
    rethrow(err);
end
end

function run_check_body(label,S,ZL,ZR,opts,seed,vs_sosvar,rowsL,rowsR)
prog0 = piesos_program(S.all);
[prog,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog0,ZL,ZR,opts);
if isempty(ZR),     ZR = ZL;    rowsR = rowsL;      end
rng(seed);
d0map = random_decvars(Pmat);
[Vt,Vl] = check_V(V,Pmat,ZopL,ZopR,ZL,ZR,S.X,d0map,rowsL,rowsR);
assert(abs(Vt-Vl)<=1e-9*(1+abs(Vt)),'%s: V(x)=%.10g differs from <ZopL*zL,Pmat*ZopR*zR>=%.10g.',label,Vl,Vt);
assert(abs(Vt)>1e-6,'%s: the value %g of V is too small to test.',label,Vt);
if vs_sosvar
    osv = rmfield(opts,'pos');
    [ps,Vs,Ps,Zs] = piesos_sosvar(prog0,ZL,osv);
    assert(isequal(prog,ps),'%s: program differs from that of piesos_sosvar.',label);
    assert(isequal(size(Pmat),size(Ps)),'%s: size of Pmat differs from that of piesos_sosvar.',label);
    for k=1:numel(Ps)
        assert(isequal(Pmat{k}.C,Ps{k}.C) && isequal(Pmat{k}.dvarname,Ps{k}.dvarname),'%s: Pmat{%d} differs from that of piesos_sosvar.',label,k);
    end
    assert(isequal(ZopL,Zs) && isequal(ZopR,Zs),'%s: Zop differs from that of piesos_sosvar.',label);
    [~,Vls] = check_V(Vs,Ps,Zs,Zs,ZL,ZL,S.X,d0map,rowsL,rowsL);
    assert(abs(Vls-Vt)<=1e-9*(1+abs(Vt)),'%s: V of piesos_sosvar (%.10g) differs from %.10g.',label,Vls,Vt);
end
end

function Y = apply_side(Zop,Z,X)
% Zop applied to the entries of a monomial basis Z (identity coefficients):
% block i acts on the monomial of row i of Z.C, the product of the states in X
% listed in the row of Z.Z that has a block in row i
Zf = full(Z.Z);
nr = size(Z.C.ops,1);
Xc = cell(1,nr);
for i=1:nr
    k = find(~cellfun(@isempty,Z.C.ops(i,:)));
    assert(isscalar(k),'Row %d of the basis does not have a single monomial.',i);
    idc = Zf(k,Zf(k,:)~=0);
    if isempty(idc)
        Xc{i} = {1};
    else
        Xc{i} = cellfun(@(nm) X.(nm),reshape(Z.varname(idc),1,[]),'UniformOutput',false);
    end
end
Y = apply_op(fix_empty_params_zop(Zop),Xc);
end

function Y = eval_rows(Zop,rows)
% Zop applied to the entries of a vector of distributed polynomials, given by
% rows{i}: the terms of entry i, each a 1 x d cell with the factors of the term
% (polynomials, the outputs of the coefficients of the states, as computed with
% apply_sopvar), to which the tensor product, the basis operator of row i, is
% applied as to the factors of a monomial; the entry is the sum over the terms.
% This does not use the composition of the operators in piesos_quadvar
Zop = fix_empty_params_zop(Zop);
Yc = cell(numel(rows),1);
for i=1:numel(rows)
    Zi = tensopmat({Zop.ops{i,i}});
    Yi = polynomial(zeros(size(Zop.ops{i,i},1),1));      % no terms: zero entry
    for t=1:numel(rows{i})
        Yi = Yi + apply_op(Zi,{rows{i}{t}});
    end
    Yc{i} = Yi;
end
Y = vertcat(Yc{:});
end

function [vn,dm] = used_vars(ZL,ZR)
% Variables on which some monomial in ZL or ZR depends, and their domains
vn = cell(1,0);     dm = zeros(0,2);
for Z = {ZL,ZR}
    Zc = Z{1}.C;
    for k=1:numel(Zc.vars.in)
        if any(Zc.depmat.in(:,k)>0) && ~any(strcmp(vn,Zc.vars.in{k}))
            vn{end+1} = Zc.vars.in{k};              %#ok<AGROW>
            dm = [dm;Zc.dom.in(k,:)];               %#ok<AGROW>
        end
    end
end
end

function Y = eval_V_generic(V,vals,dep_state,vnames)
% Evaluate the functional V (polyopvar of intop kernels) at the states
% vals.(name). The dummy variables of s_k are numbered over the factors of the
% monomial that depend on s_k, and named s_k_j (as in apply_op) if the monomial
% has more than 1 factor. dep_state(i,k): state i (order of V.varname) depends on s_k.
Y = 0;
nZ = size(V.Z,1);
Zf = full(V.Z);
for k=1:nZ
    op = V.C.ops{1,k};
    if isempty(op), continue, end
    idcs = Zf(k,Zf(k,:)~=0);
    d = numel(idcs);
    if d==0
        if isa(op,'intop')
            Y = Y + apply_op(op,polynomial(1));
        else
            Y = Y + apply_op(op,{1});
        end
    elseif isa(op,'intop')
        slot = zeros(1,numel(vnames));
        v = 1;
        for j=1:d
            xj = vals.(V.varname{idcs(j)});
            for mm=1:numel(vnames)
                if dep_state(idcs(j),mm)
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
        for j=1:d
            xv{j} = vals.(V.varname{idcs(j)});
        end
        Y = Y + apply_op(op,xv);
    end
end
end

%% ---- workaround for apply_sopvar (as in Test_piesos_sosvar.m) ---------
function Zop2 = fix_empty_params_zop(Zop)
% Expand params{k}=[] (excluded combinations) to zero sparse matrices
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


%% ---- cases ------------------------------------------------------------
function case_equals_sosvar()
% Equal bases and opts.pos=true: program, Pmat, Zop and V of piesos_sosvar (for
% scalar states, for which it is verified), with ZR omitted, empty, or equal to ZL
S = make_states();
cs = {};
cs{end+1} = {'x1, sep',S.x1,struct('degs',struct('in',2,'out',1),'sep',true)};
cs{end+1} = {'[x1;x1 o x1]',[S.x1;otimes(S.x1,S.x1)],struct('degs',struct('in',0,'out',1))};
cs{end+1} = {'[x1;x4;x1 o x4], cross term',[S.x1;S.x4;otimes(S.x1,S.x4)],struct('degs',struct('in',0,'out',0))};
cs{end+1} = {'[x1;x2], mixed dependence',[S.x1;S.x2],struct('degs',struct('in',0,'out',1))};
cs{end+1} = {'[x1;x3], disjoint variables',[S.x1;S.x3],struct('degs',struct('in',1,'out',0))};
cs{end+1} = {'x2, integrals in both variables',S.x2,struct('degs',struct('in',1,'out',0),'exclude',[true;false;true])};
cs{end+1} = {'x1, psatz 1',S.x1,struct('degs',struct('in',1,'out',2),'psatz',1)};
cs{end+1} = {'[x1;x2], psatz [0 1], exclude',[S.x1;S.x2],struct('degs',struct('in',0,'out',[1 2]),'psatz',[0 1],'exclude',[true;false;false])};
for c=1:numel(cs)
    [label,Z,opts] = cs{c}{:};
    opts.pos = true;
    run_check(label,S,Z,[],opts,c,true);
    prog0 = piesos_program(S.all);
    A = cell(1,5);      B = cell(1,5);
    [A{:}] = piesos_quadvar(prog0,Z,[],opts);
    [B{:}] = piesos_quadvar(prog0,Z,Z,opts);
    assert_same_output(A,B,[label,': ZR=ZL']);
end
end

function case_general_equal_bases()
% Equal bases, general or symmetric Pmat: all pairs of monomials (no
% shortcut for a general Pmat), V against the direct evaluation
S = make_states();
cs = {};
cs{end+1} = {'x1',S.x1,struct('degs',struct('in',1,'out',1))};
cs{end+1} = {'[x1;x1 o x1]',[S.x1;otimes(S.x1,S.x1)],struct('degs',struct('in',0,'out',1))};
cs{end+1} = {'[x1;x4;x1 o x4]',[S.x1;S.x4;otimes(S.x1,S.x4)],struct('degs',struct('in',0,'out',0))};
cs{end+1} = {'[x1;x2]',[S.x1;S.x2],struct('degs',struct('in',0,'out',0))};
cs{end+1} = {'[x1;x3]',[S.x1;S.x3],struct('degs',struct('in',0,'out',1))};
for c=1:numel(cs)
    [label,Z,opts] = cs{c}{:};
    run_check([label,' general'],S,Z,[],opts,10+c);
    run_check([label,' sym'],S,Z,[],setfield(opts,'sym',true),20+c);       %#ok<SFLD>
end
end

function case_different_bases()
% Different bases on the left and right: other monomials, states (different
% lists, dependence on different variables, constants), and numbers of them;
% general Pmat, and symmetric or SOS if the blocks match
S = make_states();
cn = constant2polyopvar(1);
cs = {};
dg = @(i,o) struct('degs',struct('in',i,'out',o));
cs{end+1} = {'ZL=[x1;x1 o x1], ZR=x1',[S.x1;otimes(S.x1,S.x1)],S.x1,dg(1,0)};
cs{end+1} = {'ZL=x1, ZR=[x1;x1 o x1]',S.x1,[S.x1;otimes(S.x1,S.x1)],dg(1,1)};
cs{end+1} = {'ZL=x1, ZR=x4',S.x1,S.x4,dg(1,1)};
cs{end+1} = {'ZL=x1, ZR=x2, nested dependence',S.x1,S.x2,dg(0,1)};
cs{end+1} = {'ZL=x2, ZR=x1, nested dependence',S.x2,S.x1,dg(0,0)};
cs{end+1} = {'ZL=[x1;x2], ZR=x3, disjoint',[S.x1;S.x2],S.x3,dg(0,1)};
cs{end+1} = {'ZL=x1, ZR=[x1;x2;x3]',S.x1,[S.x1;S.x2;S.x3],dg(0,0)};
cs{end+1} = {'ZL=[x1;x4], ZR=[x4;x1;x1 o x4], one state list in a different order',[S.x1;S.x4],[S.x4;S.x1;otimes(S.x1,S.x4)],dg(0,0)};
cs{end+1} = {'ZL=[1;x1], ZR=x1, constant monomial',[cn;S.x1],S.x1,dg(0,1)};
cs{end+1} = {'ZL=[1;x1], ZR=[1;x2], constants in both',[cn;S.x1],[cn;S.x2],dg(0,0)};
cs{end+1} = {'ZL=x2, ZR=x5, integrals in both variables',S.x2,S.x5,struct('degs',struct('in',1,'out',0),'exclude',[true;false;true])};
for c=1:numel(cs)
    [label,ZL,ZR,opts] = cs{c}{:};
    run_check(label,S,ZL,ZR,opts,30+c);
end
% symmetric and SOS Pmat for bases of the same blocks, ZL~=ZR (x1, x4 have the
% same dimension and dependence), and for the permuted [x1;x4]
for md = {'sym','pos'}
    opts = dg(1,1);     opts.(md{1}) = true;
    run_check(['ZL=x1, ZR=x4, ',md{1}],S,S.x1,S.x4,opts,50);
    run_check(['ZL=[x1;x4], ZR=[x4;x1], ',md{1}],S,[S.x1;S.x4],[S.x4;S.x1],opts,51);
    run_check(['ZL=[x1;x1 o x1], ZR=[x4;x4 o x4], ',md{1}],S,[S.x1;otimes(S.x1,S.x1)],[S.x4;otimes(S.x4,S.x4)],dg_with(opts,0,1),52);
    opts2 = struct('degs',struct('in',1,'out',0),'exclude',[true;false;true]);     opts2.(md{1}) = true;
    run_check(['ZL=x2, ZR=x5, integrals in both variables, ',md{1}],S,S.x2,S.x5,opts2,53);
end
end

function opts = dg_with(opts,i,o)
opts.degs = struct('in',i,'out',o);
end

function case_operators_differ()
% Same basis on both sides, symmetric or SOS Pmat, but ZopL~=ZopR (lower
% integral excluded on the left, upper on the right): terms (i,j) and (j,i)
% then differ, so all pairs must be used
S = make_states();
for md = {'sym','pos'}
    opts = struct('degs',struct('in',1,'out',1),'L',struct('exclude',[false;true;false]),...
        'R',struct('exclude',[false;false;true]));
    opts.(md{1}) = true;
    for Z = {S.x1,[S.x1;S.x4],[S.x1;otimes(S.x1,S.x1)]}
        run_check(['ZopL~=ZopR, ',md{1}],S,Z{1},[],opts,60);
        prog0 = piesos_program(S.all);
        [~,~,~,ZL,ZR] = piesos_quadvar(prog0,Z{1},[],opts);
        assert(~isequal(ZL,ZR),'ZopL and ZopR must differ.');
    end
end
end

function case_vector_valued()
% Vector-valued states: the kernels of the monomials are on kron(x,y), so the
% order of the factors matters for the components, for symmetric Pmat too
S = make_states();
dg = @(i,o) struct('degs',struct('in',i,'out',o));
cs = {};
cs{end+1} = {'y, pos',S.y,[],setfield(dg(1,0),'pos',true)};                %#ok<SFLD>
cs{end+1} = {'y, sym',S.y,[],setfield(dg(0,1),'sym',true)};                %#ok<SFLD>
cs{end+1} = {'y, general',S.y,[],dg(1,0)};
cs{end+1} = {'[x1;y], pos',[S.x1;S.y],[],setfield(dg(0,0),'pos',true)};    %#ok<SFLD>
cs{end+1} = {'[x1;y], general',[S.x1;S.y],[],dg(0,1)};
cs{end+1} = {'ZL=y, ZR=x1',S.y,S.x1,dg(1,0)};
cs{end+1} = {'ZL=[x1;y], ZR=y',[S.x1;S.y],S.y,dg(0,0)};
cs{end+1} = {'ZL=y, ZR=[x1;y], general',S.y,[S.x1;S.y],dg(0,0)};
% monomials with vector-valued factors: y o y, of size 4 for factors of size 2
cs{end+1} = {'[y;y o y], general',[S.y;otimes(S.y,S.y)],[],setfield(dg(0,0),'exclude',[true;true;false])};     %#ok<SFLD>
cs{end+1} = {'ZL=y o y, ZR=[y;y o y], general',otimes(S.y,S.y),[S.y;otimes(S.y,S.y)],setfield(dg(0,0),'exclude',[true;true;false])};     %#ok<SFLD>
for c=1:numel(cs)
    [label,ZL,ZR,opts] = cs{c}{:};
    run_check(label,S,ZL,ZR,opts,70+c);
end
end

function case_basis_forms()
% ZR omitted, empty, equal to ZL, or a separately declared equal basis, gives
% the same output; also for non-symmetric Pmat
S = make_states();
Z = [S.x1;S.x2];
Zcopy = [S.x1;S.x2];
prog0 = piesos_program(S.all);
for opts = {struct('pos',true,'degs',0),struct('sym',true,'degs',0),struct('degs',0)}
    A = cell(1,5);      B = A;      C = A;
    [A{:}] = piesos_quadvar(prog0,Z,[],opts{1});
    [B{:}] = piesos_quadvar(prog0,Z,Z,opts{1});
    [C{:}] = piesos_quadvar(prog0,Z,Zcopy,opts{1});
    assert_same_output(A,B,'ZR=ZL');
    assert_same_output(A,C,'ZR=copy of ZL');
end
% ZR and opts omitted
A = cell(1,5);      B = A;      C = A;
[A{:}] = piesos_quadvar(prog0,Z);
[B{:}] = piesos_quadvar(prog0,Z,[]);
[C{:}] = piesos_quadvar(prog0,Z,[],struct());
assert_same_output(A,B,'opts omitted');
assert_same_output(A,C,'ZR and opts omitted');
% bases with different states: V is in the states of both
[~,V,Pmat] = piesos_quadvar(prog0,S.x1,S.x3);
assert(isequal(size(Pmat),[1,1]) && isequal(reshape(V.varname,1,[]),{'x1','x3'}),'V must be in the states of both ZL and ZR.');
end

function case_degs_forms()
% The degrees of the left and right basis in opts.degs, or in opts.L.degs and
% opts.R.degs (instead of opts.degs.L and opts.degs.R): the same output
S = make_states();
Z = [S.x1;S.x4];
prog0 = piesos_program(S.all);
dL = struct('in',1,'out',0);        dR = struct('in',{{0,1}},'out',0);
ref = cell(1,5);
[ref{:}] = piesos_quadvar(prog0,Z,[],struct('degs',struct('L',dL,'R',dR)));
assert(~isequal(Zdims(ref{4}),Zdims(ref{5})),'The test needs different degrees for ZopL, ZopR.');
forms = {struct('L',struct('degs',dL),'R',struct('degs',dR)),...
         struct('degs',struct('L',dL),'R',struct('degs',dR)),...
         struct('degs',struct('R',dR),'L',struct('degs',dL)),...
         struct('degs',struct('L',dL,'R',dR),'exclude',[false;false;false]),...
         struct('L',struct('degs',dL),'R',struct('degs',dR),'exclude',[false;false;false])};
for k=1:numel(forms)
    out = cell(1,5);
    [out{:}] = piesos_quadvar(prog0,Z,[],forms{k});
    assert_same_output(ref,out,sprintf('degs form %d',k));
end
% the same degrees for both: opts.degs, a single field, or opts.L.degs only
A = cell(1,5);
[A{:}] = piesos_quadvar(prog0,Z,[],struct('degs',dL));
for opts = {struct('degs',struct('L',dL)),struct('degs',struct('R',dL)),struct('L',struct('degs',dL)),...
            struct('R',struct('degs',dL)),struct('degs',dL,'L',struct('degs',dL))}
    out = cell(1,5);
    [out{:}] = piesos_quadvar(prog0,Z,[],opts{1});
    assert_same_output(A,out,'the same degrees for both');
end
% a general opts.degs with a specific opts.L.degs or opts.R.degs
out = cell(1,5);
[out{:}] = piesos_quadvar(prog0,Z,[],struct('degs',dR,'L',struct('degs',dL)));
assert_same_output(ref,out,'opts.degs and opts.L.degs');
out = cell(1,5);
[out{:}] = piesos_quadvar(prog0,Z,[],struct('degs',dL,'R',struct('degs',dR)));
assert_same_output(ref,out,'opts.degs and opts.R.degs');
% value check with different degrees
S2 = make_states();
run_check('different degrees in opts.L.degs, opts.R.degs',S2,[S2.x1;S2.x4],[],struct('L',struct('degs',dL),'R',struct('degs',dR)),80);
end

function case_operator_options()
% The other options of the operators: opts.exclude etc. apply to both if opts.L
% and opts.R are absent, or only specify degrees; opts.L, opts.R override
S = make_states();
prog0 = piesos_program(S.all);
rows = @(opts) deal_rows(prog0,S.x1,opts);
ex = [true;false;false];
dg = struct('in',1,'out',0);
assert(isequal(rows(struct('degs',dg)),[5,5]),'Default rows.');
assert(isequal(rows(struct('degs',dg,'exclude',ex)),[4,4]),'Top-level exclude applies to both.');
assert(isequal(rows(struct('exclude',ex,'L',struct('degs',dg),'R',struct('degs',dg))),[4,4]),...
    'Top-level exclude applies to both if opts.L, opts.R only specify degrees.');
assert(isequal(rows(struct('exclude',ex,'L',struct('degs',dg))),[4,4]),...
    'Top-level exclude applies to both if only opts.L specifies degrees.');
assert(isequal(rows(struct('L',struct('degs',dg,'exclude',ex))),[4,4]),'Only opts.L with an option: used for both.');
assert(isequal(rows(struct('L',struct('degs',dg,'exclude',ex),'R',struct('degs',dg))),[4,5]),...
    'opts.L with an option, opts.R only degrees: the defaults for ZopR.');
assert(isequal(rows(struct('L',struct('degs',dg,'exclude',ex),'R',struct('degs',dg,'exclude',[false;true;false]))),[4,3]),...
    'opts.L and opts.R.');
assert(isequal(rows(struct('exclude',ex,'L',struct('degs',dg,'sep',true))),[3,3]),...
    'opts.L with an option (sep: 3 rows) is used for both, and the top-level exclude is ignored.');
end

function r = deal_rows(prog0,Z,opts)
% Row counts [rows of ZopL, rows of ZopR] of the operators for the basis Z
[~,~,~,ZL,ZR] = piesos_quadvar(prog0,Z,[],opts);
r = [size(ZL,1),size(ZR,1)];
end

function case_wiring()
% Pmat, ZopL, ZopR and the program must be those of piesos_quadopvar for the
% function spaces of the monomials written out by hand (the input spaces of
% ZL.C, ZR.C)
S = make_states();
prog0 = piesos_program(S.all);
cs = {{[S.x1;otimes(S.x1,S.x1)],[],struct('degs',struct('in',0,'out',1),'pos',true)},...
      {[S.x1;S.x2],S.x3,struct('degs',struct('in',1,'out',0))},...
      {S.y,S.x1,struct('degs',struct('in',1,'out',0))},...
      {[S.x1;S.x2],[],struct('degs',struct('L',struct('in',1,'out',0),'R',struct('in',0,'out',0)))}};
for c=1:numel(cs)
    [ZL,ZR,opts] = cs{c}{:};
    A = cell(1,5);
    [A{:}] = piesos_quadvar(prog0,ZL,ZR,opts);
    if isempty(ZR), ZR = ZL; end
    % by hand: dimensions, variables, domains, dependence of the monomials
    mL = double(ZL.C.depmat.in);        mR = double(ZR.C.depmat.in);
    dims = struct('L',ZL.C.dims.in,'R',ZR.C.dims.in);
    vars = struct('L',{ZL.C.vars.in},'R',{ZR.C.vars.in});
    dom = struct('L',ZL.C.dom.in,'R',ZR.C.dom.in);
    depmat = struct('L',mL,'R',mR);
    op2 = rmfield(opts,'degs');
    [p2,P2,ZL2,ZR2] = piesos_quadopvar(prog0,dims,vars,dom,depmat,opts.degs,op2);
    assert(isequal(A{1},p2),'Case %d: program differs from piesos_quadopvar.',c);
    assert(isequal(size(A{3}),size(P2)),'Case %d: Pmat has another size than that of piesos_quadopvar.',c);
    for k=1:numel(P2)
        assert(isequal(A{3}{k}.C,P2{k}.C) && isequal(A{3}{k}.dvarname,P2{k}.dvarname),'Case %d: Pmat{%d} differs.',c,k);
    end
    assert(isequal(A{4},ZL2) && isequal(A{5},ZR2),'Case %d: ZopL or ZopR differs from piesos_quadopvar.',c);
    assert(isequal(size(A{3}),[numel(dims.L),numel(dims.R)]),'Case %d: wrong number of blocks of Pmat.',c);
end
end

function case_coefficients()
% Coefficients ZL.C, ZR.C other than the identity: V is the quadratic form of
% the entries of ZL, ZR, with the operators acting on the output of their
% coefficients, e.g. on T*x1, or on the tensor product (T*x1) o (T*x1). The
% ground truth builds the entries from the outputs f=T*x (apply_sopvar) of the
% coefficients, applies the basis operator of the row to the tensor product of
% these (eval_rows), and has Pmat indexed by the rows: it does not use the
% composition of the operators in piesos_quadvar
S = make_states();
X = S.X;
rng(7);
vars = struct('in',{{'s1'}},'out',{{'s1'}});
dom = struct('in',S.dom(1,:),'out',S.dom(1,:));
T = rand_sopvar([1,1],vars,dom,1);                              % PI operators
T2 = rand_sopvar([1,1],vars,dom,1);
Tm = rand_sopvar([2,2],vars,dom,1);                             % on the vector state y
fm = apply_sopvar(Tm,X.y);
f1 = apply_sopvar(T,X.x1);              f2 = apply_sopvar(T2,X.x4);
g2 = apply_sopvar(T2,X.x1);
y1 = X.y(1);    y2 = X.y(2);
dg = @(i,o) struct('degs',struct('in',i,'out',o));
cs = {};
% label, ZL, ZR (empty: same), opts, factors in the rows of ZL (rowsL), of ZR
cs{end+1} = {'2*x1, scaling',2*S.x1,[],dg(1,0),{{{2*X.x1}}}};
cs{end+1} = {'T*x1, PI operator',T*S.x1,[],dg(1,0),{{{f1}}}};
cs{end+1} = {'T*x1, PI operator, sym',T*S.x1,[],setfield(dg(0,1),'sym',true),{{{f1}}}};       %#ok<SFLD>
cs{end+1} = {'dmonomials(T*x1,[1;2]), tensor product (T o T)*(x1 o x1)',dmonomials(T*S.x1,[1;2]),[],dg(0,1),{{{f1}},{{f1,f1}}}};
cs{end+1} = {'[x4;T*x1], identity and operator',[S.x4;T*S.x1],[],dg(0,1),{{{X.x4}},{{f1}}}};
cs{end+1} = {'T*x1+T2*x4, row with two monomials',T*S.x1+T2*S.x4,[],dg(1,0),{{{f1},{f2}}}};
f12 = {{f1,f1},{f1,f2},{f2,f1},{f2,f2}};
% (the cost of the products of sums is that of the pairs of terms: upper integral only)
cs{end+1} = {'dmonomials(T*x1+T2*x4,2), product of sums',dmonomials(T*S.x1+T2*S.x4,2),[],setfield(dg(0,0),'exclude',[true;true;false]),{f12}};   %#ok<SFLD>
cs{end+1} = {'otimes(T*x1,T*x1), tensor product of the operator outputs',[T*S.x1;otimes(T*S.x1,T*S.x1)],[],setfield(dg(0,0),'exclude',[true;false;false]),{{{f1}},{{f1,f1}}}};   %#ok<SFLD>
cs{end+1} = {'[1 2;0 3]*y, one vector-valued entry',[1 2;0 3]*S.y,[],dg(1,0),{{{[1 2;0 3]*X.y}}}};
cs{end+1} = {'[[1 2]*y;[0 3]*y], two scalar entries of one monomial',[[1 2]*S.y;[0 3]*S.y],[],dg(1,0),{{{[1 2]*X.y}},{{[0 3]*X.y}}}};
% scalar entries of the vector state y and their tensor products, by selection of its rows (dmonomials
% has the variables y, not its entries: dmonomials(y,[1 0;0 1;2 0;1 1;0 2]) is not allowed)
sy = {S.y(1,:),S.y(2,:)};
cs{end+1} = {'[y(1);y(2);y(1)o y(1);y(1)o y(2);y(2)o y(2)], selections, products',[sy{1};sy{2};otimes(sy{1},sy{1},2);otimes(sy{1},sy{2},2);otimes(sy{2},sy{2},2)],[],setfield(dg(0,0),'exclude',[true;true;false]),...   %#ok<SFLD>
             {{{y1}},{{y2}},{{y1,y1}},{{y1,y2}},{{y2,y2}}}};
% monomials of the vector-valued variable y, of degree 1 and 2: y and y o y of size 4 (n^d for n=2)
cs{end+1} = {'dmonomials(y,[1;2]), vector variable',dmonomials(S.y,[1;2]),[],setfield(dg(0,0),'exclude',[true;true;false]),{{{X.y}},{{X.y,X.y}}}};   %#ok<SFLD>
cs{end+1} = {'dmonomials(Tm*y,[1;2]), operator on a vector variable',dmonomials(Tm*S.y,[1;2]),[],setfield(dg(0,0),'exclude',[true;true;false]),{{{fm}},{{fm,fm}}}};   %#ok<SFLD>
% vector-valued factors: (Tm*y) o (Tm*y) of size 4 for Tm*y of size 2
cs{end+1} = {'[Tm*y;otimes(Tm*y,Tm*y)], vector factors',[Tm*S.y;otimes(Tm*S.y,Tm*S.y)],[],setfield(dg(0,0),'exclude',[true;true;false]),{{{fm}},{{fm,fm}}}};   %#ok<SFLD>
cs{end+1} = {'ZL=T*x1, ZR=[x4;T2*x1]',T*S.x1,[S.x4;T2*S.x1],dg(1,0),{{{f1}}},{{{X.x4}},{{g2}}}};
cs{end+1} = {'ZL=dmonomials(T*x1,[1;2]), ZR=T2*x4',dmonomials(T*S.x1,[1;2]),T2*S.x4,dg(0,1),{{{f1}},{{f1,f1}}},{{{f2}}}};
cs{end+1} = {'ZL=[x1;x1 o x1], ZR=dmonomials(T*x1,[1;2]), identity and operator',[S.x1;otimes(S.x1,S.x1)],dmonomials(T*S.x1,[1;2]),dg(0,0),...
             {{{X.x1}},{{X.x1,X.x1}}},{{{f1}},{{f1,f1}}}};
for c=1:numel(cs)
    [label,ZL,ZR,opts,rowsL] = cs{c}{1:5};
    if numel(cs{c})>=6,     rowsR = cs{c}{6};   else,   rowsR = [];     end
    tc = tic;
    try
        run_check(label,S,ZL,ZR,opts,90+c,false,rowsL,rowsR);
    catch err
        error('%s: %s',label,err.message);
    end
    fprintf('        %s (%.1f s)\n',label,toc(tc));
end
% excluded terms (empty parameters of ZopL, ZopR must not break the composition)
for ex = {[true;false;false],[false;true;false],[false;true;true]}
    op_ex = setfield(dg(0,0),'exclude',ex{1});                  %#ok<SFLD>
    run_check(['dmonomials(T*x1,[1;2]), exclude ',mat2str(ex{1}')],S,dmonomials(T*S.x1,[1;2]),[],op_ex,110,false,{{{f1}},{{f1,f1}}});
end
% positive Pmat for ZL=ZR with an operator: the same as piesos_sosvar
run_check('T*x1, pos',S,T*S.x1,[],struct('degs',struct('in',1,'out',1),'pos',true),120,true,{{{f1}}});
run_check('dmonomials(T*x1,[1;2]), pos',S,dmonomials(T*S.x1,[1;2]),[],struct('degs',struct('in',0,'out',1),'pos',true,'exclude',[true;true;false]),121,true,{{{f1}},{{f1,f1}}});
% Pmat is indexed by the rows of ZL.C, ZR.C: duplicates and unused monomials
prog0 = piesos_program(S.all);
[~,~,P] = piesos_quadvar(prog0,[S.x1;S.x1],[],dg(0,0));
assert(isequal(size(P),[2,2]),'Pmat must have a block for each row of ZL.C, ZR.C (duplicate rows).');
[~,~,P] = piesos_quadvar(prog0,dmonomials([S.x1;S.x4],(1:2)'),[],dg(0,0));
assert(isequal(size(P),[5,5]),'Pmat must have a block for each row of ZL.C, ZR.C (unused monomial).');
[~,~,P] = piesos_quadvar(prog0,dmonomials(T*S.x1,[1;2]),T2*S.x4,dg(0,0));
assert(isequal(size(P),[2,1]),'Pmat must have a block for each row of ZL.C, ZR.C (rectangular).');
% a zero entry (0*x1): no contribution to V, but a block in Pmat
[~,~,P] = piesos_quadvar(prog0,[S.x1;0*S.x1],[],dg(0,0));
assert(isequal(size(P),[2,2]),'Pmat must have a block for each row of ZL.C, ZR.C (zero entry).');
run_check('[x1;0*x1], zero entry',S,[S.x1;0*S.x1],[],dg(1,0),130,false,{{{X.x1}},{}});
% an operator on the coefficients changes V: the oracle with the identity in
% place of T gives another value, so a dropped coefficient would be detected
assert_sensitive('T*x1',S,T*S.x1,dg(1,0),{{{f1}}},{{{X.x1}}},131);
assert_sensitive('dmonomials(T*x1,[1;2])',S,dmonomials(T*S.x1,[1;2]),dg(0,1),{{{f1}},{{f1,f1}}},{{{X.x1}},{{X.x1,X.x1}}},132);
end

function assert_sensitive(label,S,ZL,opts,rowsL,rowsId,seed)
% The value of V for the factors in rowsL differs from that with the factors
% in rowsId (the states themselves, i.e. the coefficients dropped)
prog0 = piesos_program(S.all);
[~,V,Pmat,ZopL,ZopR] = piesos_quadvar(prog0,ZL,[],opts);
rng(seed);
d0map = random_decvars(Pmat);
[Vt,Vl] = check_V(V,Pmat,ZopL,ZopR,ZL,ZL,S.X,d0map,rowsL,rowsL);
Vid = check_V(V,Pmat,ZopL,ZopR,ZL,ZL,S.X,d0map,rowsId,rowsId);
assert(abs(Vt-Vl)<=1e-9*(1+abs(Vt)),'%s: V differs from the direct evaluation.',label);
assert(abs(Vt-Vid)>1e-6*(1+abs(Vt)),'%s: the test cannot detect dropped coefficients (%.6g and %.6g).',label,Vt,Vid);
end

function case_ode_state()
% A state without spatial variables (ODE state) next to a PDE state: the output
% space of its coefficients has no variables (empty dependence in
% coefficient_blocks), also with a coefficient
S = make_states();
xo = polyopvar({'xo'},1,{'s1'},S.dom(1,:),false);
S.all = [S.x1;xo];
S.X.xo = polynomial(1.3);
rng(7);
vars = struct('in',{{'s1'}},'out',{{'s1'}});
dom = struct('in',S.dom(1,:),'out',S.dom(1,:));
T = rand_sopvar([1,1],vars,dom,1);
f1 = apply_sopvar(T,S.X.x1);
dg = @(i,o) struct('degs',struct('in',i,'out',o));
run_check('[xo;x1]',S,[xo;S.x1],[],dg(0,1),140);
% a product of ODE states in the basis is not supported (no copy of a variable)
assert_err(@() piesos_quadvar(piesos_program(S.all),[xo;otimes(xo,xo);S.x1],[],dg(0,1)),'without spatial variables');
run_check('3*xo, scaling',S,3*xo,[],dg(0,0),142,false,{{{3*S.X.xo}}});
run_check('[xo;T*x1]',S,[xo;T*S.x1],[],dg(0,1),143,false,{{{S.X.xo}},{{f1}}});
end

function case_coefficient_errors()
% Coefficients that are not supported are rejected with an error, not dropped
S = make_states();
prog0 = piesos_program(S.all);
opts = struct('degs',0);
rng(7);
vars = struct('in',{{'s1'}},'out',{{'s1'}});
dom = struct('in',S.dom(1,:),'out',S.dom(1,:));
T = rand_sopvar([1,1],vars,dom,1);
% an elementwise product in the basis (not the tensor product of the variables)
assert_err(@() piesos_quadvar(prog0,[S.x1;S.x1.*S.x1],[],opts),'tensor products of the variables');
assert_err(@() piesos_quadvar(prog0,S.x1,[S.x1;S.x1.*S.x1],opts),'basis ZR must be defined in terms of tensor products');
assert_err(@() piesos_sosvar(prog0,[S.x1;S.x1.*S.x1],opts),'tensor products of the variables');
% a pointwise product nested in a tensor product (set by hand, so that the test
% does not depend on how dmonomials builds products)
Zn = [S.x1;S.x1.*S.x4];
[i2,j2] = find(cellfun(@(b) isa(b,'tensopvar') && size(b.ops,2)==2,Zn.C.ops));
assert(~isempty(i2),'The test needs a block with two factors.');
blk = Zn.C.ops{i2(1),j2(1)};    blk.ops = {blk};
Zn.C.ops{i2(1),j2(1)} = blk;    Zn.chkval = true;
assert_err(@() piesos_quadvar(prog0,Zn,[],opts),'nested factors');
assert_err(@() piesos_quadvar(prog0,S.x1,Zn,opts),'ZR with nested factors');
% a zero row
Z0 = [S.x1;S.x4];
Z0.C.ops(2,:) = cell(1,size(Z0.C.ops,2));       Z0.chkval = true;
assert_err(@() piesos_quadvar(prog0,Z0,[],opts),'Row 2 of ZL.C is zero');
assert_err(@() piesos_quadvar(prog0,S.x1,Z0,opts),'Row 2 of ZR.C is zero');
% a monomial with factors of different sizes (not constructed by dmonomials or
% otimes: set by hand in the product y o y)
Zm = [S.y;otimes(S.y,S.y)];
[i2,j2] = find(cellfun(@(b) isa(b,'tensopvar') && size(b.ops,2)==2,Zm.C.ops));
assert(~isempty(i2),'The test needs a block with two factors.');
blk = Zm.C.ops{i2(1),j2(1)};    blk.dims(2,1) = blk.dims(1,1)+1;
Zm.C.ops{i2(1),j2(1)} = blk;    Zm.chkval = true;
assert_err(@() piesos_quadvar(prog0,Zm,[],opts),'must have the same size');
assert_err(@() piesos_quadvar(prog0,S.x1,Zm,opts),'must have the same size');
% a block that is not a tensor-PI operator
Zb = T*S.x1;
Zb.C.ops{1,1} = T;                              Zb.chkval = true;
assert_err(@() piesos_quadvar(prog0,Zb,[],opts),'must be tensor-PI operators');
% a row with blocks acting on functions of different variables
Zs = [S.x1;S.x2];
Zs.C.ops{1,2} = Zs.C.ops{2,2};                  Zs.chkval = true;
assert_err(@() piesos_quadvar(prog0,Zs,[],opts),'must match that of output');
end

function case_errors()
S = make_states();
prog0 = piesos_program(S.all);
assert_err(@() piesos_quadvar(prog0,5),"'polyopvar'");
assert_err(@() piesos_quadvar(prog0,S.x1,5),"'polyopvar'");
assert_err(@() piesos_quadvar(prog0,5,S.x1),"'polyopvar'");
assert_err(@() piesos_quadvar(prog0,S.x1,[],5),"'struct'");
assert_err(@() piesos_quadvar(prog0,S.x1,S.x1,5),"'struct'");
% a monomial of states depending on different variables (in ZL or in ZR)
Zbad = [S.x1;S.x2;otimes(S.x1,S.x2)];
assert_err(@() piesos_quadvar(prog0,Zbad,[],struct('sep',true)),'different variables');
assert_err(@() piesos_quadvar(prog0,Zbad,S.x1,struct('sep',true)),'different variables');
assert_err(@() piesos_quadvar(prog0,S.x1,Zbad,struct('sep',true)),'different variables');
% degrees specified twice for the same basis
dg = struct('in',1,'out',0);
assert_err(@() piesos_quadvar(prog0,S.x1,[],struct('degs',struct('L',dg),'L',struct('degs',dg))),'cannot be specified both');
assert_err(@() piesos_quadvar(prog0,S.x1,[],struct('degs',struct('R',dg),'R',struct('degs',dg))),'cannot be specified both');
piesos_quadvar(prog0,S.x1,[],struct('degs',struct('R',dg),'L',struct('degs',dg)));     % not twice
% errors of piesos_quadopvar: blocks of different dimensions, psatz without pos, options
assert_err(@() piesos_quadvar(prog0,S.x1,S.y,struct('pos',true)),'Output dimensions');
assert_err(@() piesos_quadvar(prog0,S.x1,[S.x1;S.x4],struct('sym',true)),'Output dimensions');
assert_err(@() piesos_quadvar(prog0,S.x1,[],struct('L',5)),'must be specified as');
assert_err(@() piesos_quadvar(prog0,S.x1,[],struct('degs',struct('in',1))),"fields 'in' and 'out'");
assert_err(@() piesos_quadvar(prog0,S.x1,[],struct('L',struct('degs',struct('in',1)))),"fields 'in' and 'out'");
assert_err(@() piesos_quadvar(prog0,S.x1,[],struct('pos',true,'degs',struct('L',dg,'R',struct('in',1,'out',1)))),'Degrees of the primary');
% states with the same name but another size, or a variable with another domain
yy = polyopvar('x1',2,polynomial({'s1'}),S.dom(1,:));
assert_err(@() piesos_quadvar(prog0,S.x1,yy),'must match');
xx = polyopvar('x1',polynomial({'s1'}),[0 1]);
assert_err(@() piesos_quadvar(prog0,S.x1,xx),'same domain');
end
