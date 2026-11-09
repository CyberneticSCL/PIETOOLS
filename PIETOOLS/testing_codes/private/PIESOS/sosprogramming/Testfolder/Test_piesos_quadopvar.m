function Test_piesos_quadopvar()
% TEST_PIESOS_QUADOPVAR tests, using assertions, piesos_quadopvar,
%   [prog,Pmat,ZopL,ZopR] = ...
%               piesos_quadopvar(prog,dims,vars,dom,depmat,degs,opts),
% which declares Pop = ZopL'*Pmat*ZopR, with ZopL and ZopR bases of
% tensor-PI operators (kernels in the dummy variables only) on a left
% (L/out) and a right (R/in) function space, and Pmat a p x q cell of
% matrices in the primary variables: general, symmetric (opts.sym) or SOS
% (opts.pos).
%
% Ground truths are independent of the function under test:
%  - equal left and right inputs and opts.pos=true: the output must be
%    that of piesos_poslpivar (verified by Test_piesos_poslpivar): the
%    whole program, Pmat, and ZopL=ZopR=Zop;
%  - ZopL, ZopR: the number of rows is counted by hand from the monomial
%    bases (a 3-PI direction with dummy degree r gives 1 multiplier row
%    and (r+1) rows for each of the lower and upper integral; counts
%    multiply over the variables of a state), every kernel has ZL=0 and
%    the dummy monomials (0:r)' in ZR, and both map to the same space of
%    the variables on which some state of either depends. Left and right
%    differ in the cases, so that a swap of left and right is detected;
%  - Pmat: for each psatz term, a random Gram matrix Qk (general,
%    symmetric or PSD, as the mode) is assigned to the decision variables,
%    through the names returned by a direct call of sosquadvar on a fresh
%    program with monomial vectors ZsL, ZsR built by hand. Pmat(s),
%    evaluated at random points s, is compared with
%        sum_k g_k(s) * (I x ZsL(s))' * Qk * (I x ZsR(s)),
%    and the declared decision variables with those of the direct call;
%  - the forms of the inputs (arrays, structs with fields L/R or out/in, a
%    single field for both), and the errors for inappropriate inputs.
% No programs are solved.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - Test_piesos_quadopvar
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
% - piesos_poslpivar wraps piesos_quadopvar, so case_equals_poslpivar only
%   checks the forwarding of the inputs; the ground truths of the other
%   cases are independent;
% - Random Gram matrices and evaluation points depend on rng(k) at the
%   start of each case.

tests ={@case_equals_poslpivar, @case_input_forms, @case_degs_opts_forms, ...
         @case_structure_1D, @case_structure_blocks, @case_structure_spaces, ...
         @case_swap, @case_modes_same_dims, @case_psatz_LR, ...
         @case_sequential, @case_no_variables, @case_errors};
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
    error('Test_piesos_quadopvar: %d of %d cases failed.',nfail,numel(tests));
end
fprintf('Test_piesos_quadopvar: all %d cases passed.\n',numel(tests));
end


%% ---- helpers ----------------------------------------------------------
function prog = new_prog(varnames)
% Fresh program in the given spatial variables
N = numel(varnames);
x1 = polyopvar('x1',polynomial(varnames(:)),[zeros(N,1),(1:N)']);
prog = piesos_program(x1);
end

function m = nrows_3pi(r,excl)
% Rows of the basis in one 3-PI direction for dummy degree r; excl = logical
% [mult,lower,upper] of excluded terms. The multiplier has a single row.
m = ~excl(1)*1 + ~excl(2)*(r+1) + ~excl(3)*(r+1);
end

function m = Zdims(Zop)
% Row counts of the blocks of Zop
m = Zop.dims.out(:);
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
% Polynomial monomial vector in the order of zvec_num, built from the
% mixed-radix expansion of the index, independently of the code under test
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
% Evaluate the polynomial matrix P at the point vals, from the fields of the
% 'polynomial' object directly
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

function assert_err(f,needle)
try
    f();
catch err
    assert(contains(err.message,needle),'Wrong error message: "%s" (expected "%s").',err.message,needle);
    return
end
error('Expected an error containing "%s", none thrown.',needle);
end

function g = no_multiplier()
g = {@(s) 1};
end

function g = psatz_terms(dom,ps,lin)
% Multipliers prod_k (s_k-a_k)^jk*(b_k-s_k)^(qk-jk) as functions of s, for the
% terms with linear index lin (first variable fastest) in the array of size
% ps+1; ps(k)=order in variable k
g = cell(1,numel(lin));
for t=1:numel(lin)
    jj = cell(1,numel(ps)+1);
    [jj{:}] = ind2sub([ps(:)'+1,1],lin(t));
    j = cell2mat(jj(1:numel(ps)))-1;
    g{t} = @(s) prod((s(:)'-dom(:,1)').^j.*(dom(:,2)'-s(:)').^(ps-j));
end
end

function E = cut_E(E,ps)
% Exponents cut by the psatz multipliers: in the kth used variable, those
% above max(E)-ceil(ps(k)/2) are removed (ps(k)=order of the multiplier)
for i=1:numel(E)
    for k=1:numel(E{i})
        E{i}{k} = E{i}{k}(E{i}{k}<=max(E{i}{k})-ceil(ps(k)/2));
    end
end
end

function [vout,dout,used] = ref_primary(vL,dL,pL,vR,dR,pR)
% Primary variables, by explicit loops: those of the left operator, then the
% new ones of the right operator; used(k)=true if some state depends on it
vout = reshape(vL,1,[]);        dout = dL;      used = any(pL>0,1);
for k=1:numel(vR)
    idx = find(strcmp(vout,vR{k}));
    if isempty(idx)
        vout{end+1} = vR{k};                    %#ok<AGROW>
        dout = [dout;dR(k,:)];                  %#ok<AGROW>
        used(end+1) = any(pR(:,k)>0);           %#ok<AGROW>
    else
        used(idx) = used(idx) || any(pR(:,k)>0);
    end
end
used = logical(used);
end

function [m,ZR] = expected_rows(dims,dep,r,vs,vout,excl)
% Rows and dummy monomials expected for the blocks of Zop: block i has
% dims(i)*prod(rows per variable of the state) rows, and dummy monomials
% (0:r(k))' in the variables k of the state (indexed as vout). excl as in
% nrows_3pi, for a state in a single variable only
if nargin<6, excl = false(1,3); end
q = numel(dims);
m = zeros(q,1);         ZR = cell(q,1);
for i=1:q
    own = find(dep(i,:)>0);
    m(i) = dims(i)*prod(arrayfun(@(k) nrows_3pi(r(k),excl),own));
    ZR{i} = cell(1,numel(vout));
    for k=own
        ZR{i}{strcmp(vout,vs{k})} = (0:r(k))';
    end
end
end

function E = primary_E(q,o,used)
% Exponents (0:o(k))' in the used primary variables k, for each of q blocks
E = cell(1,q);
for i=1:q
    E{i} = arrayfun(@(k) (0:o(k))',find(used),'UniformOutput',false);
end
end

function e = exps_by_name(vals,names,varnames)
% Reorder vals (indexed as names) as varnames ([] if not in names)
e = cell(1,numel(varnames));
for k=1:numel(varnames)
    idx = find(strcmp(names,varnames{k}));
    if ~isempty(idx), e{k} = vals{idx}; end
end
end

function check_Zop(Zop,m_exp,ZR_exp,varnames,vout)
% The kernels have no primary monomial (ZL=0) and the dummy monomials ZR_exp{i}
% (1 x N cell of exponents of the first factor, [] for a variable on which the
% state does not depend). If vout is given, a single-factor block must have
% exactly these output variables
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

function check_sos_LR(prog,Pmat,ZopL,ZopR,spec,prog_fresh)
% Pmat = sum_k g_k*(I x ZsL)'*Qk*(I x ZsR) for random Gram matrices Qk, mode
% 'gen', 'sym' or 'pos' (PSD). spec fields: vused (used primary variables),
% dom (N x 2), EL (1 x p), ER (1 x q) cells of 1 x N cells of exponents of
% ZsL{i}, ZsR{j}, g (1 x T cell of handles g(s), in order of declaration), mode.
p = numel(spec.EL);     q = numel(spec.ER);
T = numel(spec.g);      N = numel(spec.vused);
assert(iscell(Pmat) && isequal(size(Pmat),[p,q]),'Pmat must be a %d x %d cell.',p,q);
mL = zeros(p,1);        nzL = zeros(p,1);       ZsL = cell(p,1);
for i=1:p
    mL(i) = size(ZopL.ops{i,i},1);
    nzL(i) = prod(cellfun(@numel,spec.EL{i}));
    ZsL{i} = zvec_poly(spec.EL{i},spec.vused);
end
mR = zeros(q,1);        nzR = zeros(q,1);       ZsR = cell(q,1);
for j=1:q
    mR(j) = size(ZopR.ops{j,j},1);
    nzR(j) = prod(cellfun(@numel,spec.ER{j}));
    ZsR{j} = zvec_poly(spec.ER{j},spec.vused);
end
for i=1:p
    for j=1:q
        assert(isequal(size(Pmat{i,j}),[mL(i),mR(j)]),'Pmat{%d,%d} has size %s, expected %s.',...
            i,j,mat2str(size(Pmat{i,j})),mat2str([mL(i),mR(j)]));
        assert(isa(Pmat{i,j},'dpvar'),'Pmat{%d,%d} must be a dpvar.',i,j);
    end
end
% Direct declarations in a fresh program, in the order of the terms
opt = {};
if strcmp(spec.mode,'pos')
    opt = {'pos'};
elseif strcmp(spec.mode,'sym')
    opt = {'sym'};
end
prog0 = prog_fresh;
Qn = cell(1,T);
for t=1:T
    [prog0,~,Qn{t}] = sosquadvar(prog0,ZsL,ZsR,mL,mR,opt{:});
end
assert(prog.var.num==prog_fresh.var.num+T && isequal(prog.var.idx,prog0.var.idx),...
    'Declared decision variables differ from %d direct sosquadvar call(s).',T);
vtype = 'poly';
if strcmp(spec.mode,'pos'), vtype = 'sos'; end
assert(all(strcmp(prog.var.type(prog_fresh.var.num+1:end),vtype)),'All new variables must be declared as %s.',vtype);
% Assign random Gram matrices
nQL = sum(mL.*nzL);     nQR = sum(mR.*nzR);
offL = [0;cumsum(mL.*nzL)];     offR = [0;cumsum(mR.*nzR)];
names = cell(0,1);      vals = zeros(0,1);
Qnum = cell(1,T);
for t=1:T
    switch spec.mode
        case 'pos'
            Qnum{t} = random_psd(nQL);
        case 'sym'
            A = randn(nQL);
            Qnum{t} = (A+A')/2;
        otherwise
            Qnum{t} = randn(nQL,nQR);
    end
    for i=1:p
        for j=1:q
            nm = Qn{t}{i,j};
            vl = Qnum{t}(offL(i)+1:offL(i+1),offR(j)+1:offR(j+1));
            names = [names;nm(:)];          %#ok<AGROW>
            vals = [vals;vl(:)];            %#ok<AGROW>
        end
    end
end
[names,iu] = unique(names,'stable');
vals = vals(iu);
d0map = containers.Map(names,num2cell(vals));
% Evaluate once (eval_dpvar_at is slow for large Gram matrices)
Pevals = cell(p,q);
for i=1:p
    for j=1:q
        Pevals{i,j} = eval_dpvar_at(Pmat{i,j},d0map);
    end
end
% Compare at random points in the domain
for k=1:3
    s0 = spec.dom(:,1)' + rand(1,N).*(spec.dom(:,2)-spec.dom(:,1))';
    Pfull = zeros(sum(mL),sum(mR));
    for i=1:p
        for j=1:q
            Pnum = evalpoly(Pevals{i,j},spec.vused,s0);
            Pexp = zeros(mL(i),mR(j));
            zi = zvec_num(spec.EL{i},s0);       zj = zvec_num(spec.ER{j},s0);
            for t=1:T
                Qij = Qnum{t}(offL(i)+1:offL(i+1),offR(j)+1:offR(j+1));
                Pexp = Pexp + spec.g{t}(s0)*kron(eye(mL(i)),zi)'*Qij*kron(eye(mR(j)),zj);
            end
            assert(norm(Pnum-Pexp,'fro')<=1e-8*(1+norm(Pexp,'fro')),...
                'Pmat{%d,%d}(s) differs from the expression (error %g).',i,j,norm(Pnum-Pexp,'fro'));
            r1 = sum(mL(1:i-1))+1;      c1 = sum(mR(1:j-1))+1;
            Pfull(r1:r1+mL(i)-1,c1:c1+mR(j)-1) = Pnum;
        end
    end
    if any(strcmp(spec.mode,{'sym','pos'}))
        assert(norm(Pfull-Pfull','fro')<=1e-8*(1+norm(Pfull,'fro')),'Pmat(s) is not symmetric.');
    end
    if strcmp(spec.mode,'pos')
        assert(min(eig((Pfull+Pfull')/2))>=-1e-8*(1+norm(Pfull,'fro')),...
            'Pmat(s) is not positive semidefinite at a point of the domain.');
    end
end
end

function [prog,Pmat,ZL,ZR,prog0,spec,vout,dout,used] = run_structure(cc,mode,popts)
% Declare Pmat for the configuration cc, check ZopL and ZopR against hand
% counts, and return the spec for check_sos_LR. cc: dL,dR (dims), vL,vR
% (variables), mL,mR (domains), pL,pR (depmat), rL,rR (dummy degrees, per
% variable of vL, vR), oL,oR (primary degrees, per primary variable), and
% optionally xL,xR (excluded terms [mult,lower,upper], single-variable states).
% mode: 'gen', 'sym' or 'pos'; popts: additional fields of opts
if nargin<2, mode = 'gen'; end
if nargin<3, popts = struct(); end
[vout,dout,used] = ref_primary(cc.vL,cc.mL,cc.pL,cc.vR,cc.mR,cc.pR);
prog0 = new_prog(vout);
opts = popts;
if ~strcmp(mode,'gen'), opts.(mode) = true; end
xL = false(1,3);        xR = false(1,3);
if isfield(cc,'xL')
    xL = logical(cc.xL);    xR = logical(cc.xR);
    opts.L = struct('exclude',xL(:));       opts.R = struct('exclude',xR(:));
end
[prog,Pmat,ZL,ZR] = piesos_quadopvar(prog0,struct('L',cc.dL,'R',cc.dR),...
    struct('L',{cc.vL},'R',{cc.vR}),struct('L',cc.mL,'R',cc.mR),struct('L',cc.pL,'R',cc.pR),...
    struct('L',struct('in',cc.rL,'out',cc.oL),'R',struct('in',cc.rR,'out',cc.oR)),opts);
[mLx,ZRL] = expected_rows(cc.dL,cc.pL,cc.rL,cc.vL,vout,xL);
[mRx,ZRR] = expected_rows(cc.dR,cc.pR,cc.rR,cc.vR,vout,xR);
check_Zop(ZL,mLx,ZRL,vout,vout(used));
check_Zop(ZR,mRx,ZRR,vout,vout(used));
assert(isequal(ZL.vars.out,ZR.vars.out) && isequal(ZL.dom.out,ZR.dom.out),...
    'ZopL and ZopR must map to the same space.');
spec = struct('vused',{vout(used)},'dom',dout(used,:),'EL',{primary_E(numel(cc.dL),cc.oL,used)},...
    'ER',{primary_E(numel(cc.dR),cc.oR,used)},'g',{no_multiplier()},'mode',mode);
end

function assert_same_output(A,B,label)
% Outputs A, B = {prog,Pmat,ZopL,ZopR} of two calls must be identical
assert(isequal(A{1},B{1}),'%s: programs differ.',label);
assert(isequal(size(A{2}),size(B{2})),'%s: sizes of Pmat differ.',label);
for k=1:numel(A{2})
    assert(isequal(A{2}{k}.C,B{2}{k}.C) && isequal(A{2}{k}.degmat,B{2}{k}.degmat) ...
        && isequal(A{2}{k}.varname,B{2}{k}.varname) && isequal(A{2}{k}.dvarname,B{2}{k}.dvarname) ...
        && isequal(A{2}{k}.matdim,B{2}{k}.matdim),'%s: Pmat{%d} differs.',label,k);
end
assert(isequal(A{3},B{3}),'%s: ZopL differs.',label);
assert(isequal(A{4},B{4}),'%s: ZopR differs.',label);
end

function cc = cfg(varargin)
% Configuration struct, with default dummy and primary degrees 0
cc = struct('dL',1,'dR',1,'vL',{{'s1'}},'vR',{{'s1'}},'mL',[0.5 2],'mR',[0.5 2],'pL',1,'pR',1,...
            'rL',0,'rR',0,'oL',0,'oR',0);
for k=1:2:numel(varargin)
    cc.(varargin{k}) = varargin{k+1};
end
end


%% ---- cases ------------------------------------------------------------
function case_equals_poslpivar()
% Equal left and right inputs and opts.pos=true: the output of
% piesos_poslpivar, for all formats of the inputs, psatz, exclude/sep,
% several spatial variables, mixed dependence, and no variables
rng(1);
a = 0.5;    b = 2;
d1 = [a b];     d2 = [a b;-1 0.7];      d3 = [a b;-1 0.7;1 3];
cs = {};
%       {dims, vars, dom, depmat, degs, opts}: trailing inputs may be omitted
cs{end+1} = {1,{'s1'},d1,1};
cs{end+1} = {1,{'s1'},d1,1,2};
cs{end+1} = {1,{'s1'},d1,1,struct('in',1,'out',2),struct()};
cs{end+1} = {1,{'s1'},d1,1,struct('in',{{1}},'out',{{2}}),struct()};
cs{end+1} = {1,{'s1'},d1,1,struct('in',{{{[0;1;2]}}},'out',{{{[0;2;3]}}}),struct()};
cs{end+1} = {1,{'s1'},d1,1,{{[0;1]}},struct()};
cs{end+1} = {1,{'s1'},d1,1,struct('in',2,'out',{{{[0;1;4]}}}),struct()};
cs{end+1} = {2,{'s1'}',d1,1,struct('in',1,'out',1),struct()};
cs{end+1} = {struct('in',2,'out',2),struct('in',{{'s1'}}),struct('in',[0 1]),struct('in',1),struct('in',1,'out',1),struct()};
cs{end+1} = {1,{'s1','s2'},d2,[1 1],[2 1],struct()};
cs{end+1} = {1,{'s1','s2'},d2,[1 1],struct('in',[1 2],'out',[2 0]),struct()};
cs{end+1} = {1,{'s1','s2'},d2,[1 1],struct('in',{{{[0;1],[0;2]}}},'out',{{{[0;3],[0;1;2]}}}),struct()};
cs{end+1} = {[1;2],{'s1'},d1,[1;1],struct('in',{{1,2}},'out',{{2,1}}),struct()};
cs{end+1} = {[1;1],{'s1'},d1,[1;1],struct('in',{{1,2}},'out',{{{[0;2]},{[0;1;3]}}}),struct()};
cs{end+1} = {1,{'s1','s2','s3'},d3,[1 1 1],struct('in',[1 1 0],'out',[1 0 0]),struct()};
cs{end+1} = {1,{'s1'},d1,2,struct('in',2,'out',1),struct()};
cs{end+1} = {1,{'s1'},d1,1,struct('in',2,'out',1),struct('exclude',[true;false;false])};
cs{end+1} = {1,{'s1'},d1,1,struct('in',2,'out',1),struct('sep',true)};
cs{end+1} = {1,{'s1'},d1,1,struct('in',1,'out',4),struct('psatz',2)};
cs{end+1} = {1,{'s1'},d1,1,struct('in',0,'out',3),struct('psatz',[true;false;true])};
cs{end+1} = {1,{'s1','s2'},d2,[1 1],struct('in',[1 0],'out',[2 3]),struct('psatz',[1 2])};
cs{end+1} = {1,{'s1','s2'},d2,[1 0],struct('in',1,'out',1),struct('psatz',1)};
cs{end+1} = {1,{'s2','s1'},d2([2 1],:),[1 1],struct('in',[0 0],'out',[2 1]),struct('psatz',[1 0])};
cs{end+1} = {[1;1],{'s1','s2'},d2,[1 0;1 1],struct('in',[1 1],'out',[1 2]),struct('psatz',[0 1])};
cs{end+1} = {[1;2;1],{'s1','s2'},d2,[0 0;1 0;1 1],struct('in',[1 1],'out',[1 0]),struct()};
cs{end+1} = {[1;1],{'s1','s2'},d2,[1 0;0 1],struct('in',[1 1],'out',[1 1]),struct('exclude',[false;true;false])};
cs{end+1} = {[1;1],{'s1','s2','s3'},d3,[1 0 0;1 1 0],struct('in',[1 1 1],'out',[1 1 1]),struct()};
cs{end+1} = {[2;1],{},zeros(0,2),zeros(2,0)};
for c=1:numel(cs)
    args = cs{c};
    vn = args{2};
    if isstruct(vn), vn = vn.in; end
    if isempty(vn)
        prog0 = new_prog({'s1'});
    else
        prog0 = new_prog(reshape(vn,1,[]));
    end
    [p1,P1,Z1] = piesos_poslpivar(prog0,args{:});
    % same inputs for left and right, with opts.pos
    if numel(args)<5, args{5} = 0; end
    if numel(args)<6, args{6} = struct(); end
    args{6}.pos = true;
    [p2,P2,ZL2,ZR2] = piesos_quadopvar(prog0,args{:});
    assert_same_output({p1,P1,Z1,Z1},{p2,P2,ZL2,ZR2},sprintf('Case %d',c));
end
end

function case_input_forms()
% The inputs for the left (L/out) and right (R/in) operator: a struct with
% fields L, R or the aliases out, in; an array, or a struct with only one
% field, applies to both
rng(2);
a = 0.5;    b = 2;
d2 = [a b;-1 0.7];
prog0 = new_prog({'s1','s2'});
% (a) different left and right spaces: aliases of the fields
dL = [1;2];     dR = 1;
vL = {'s1'};    vR = {'s1','s2'};
mL = [a b];     mR = d2;
pL = [1;1];     pR = [1 1];
degs = struct('L',struct('in',{{1,2}},'out',{{1,0}}),'R',struct('in',[1 2],'out',[0 1]));
ref = cell(1,4);
[ref{:}] = piesos_quadopvar(prog0,struct('L',dL,'R',dR),struct('L',{vL},'R',{vR}),...
    struct('L',mL,'R',mR),struct('L',pL,'R',pR),degs,struct());
assert(isequal(Zdims(ref{3}),[5;2*7]) && isequal(Zdims(ref{4}),5*7),'Left and right spaces are mixed up.');
vars = {struct('out',{vL},'in',{vR}),struct('L',{vL},'in',{vR}),struct('out',{vL},'R',{vR})};
for k=1:numel(vars)
    out = cell(1,4);
    [out{:}] = piesos_quadopvar(prog0,struct('out',dL,'in',dR),vars{k},struct('out',mL,'in',mR),...
        struct('L',pL,'in',pR),degs,struct());
    assert_same_output(ref,out,sprintf('(a) variant %d',k));
end
% (b) the same space for left and right
v = {2,{'s1'},[a b],1,struct('in',1,'out',2),struct()};
prog1 = new_prog({'s1'});
ref = cell(1,4);
[ref{:}] = piesos_quadopvar(prog1,v{:});
forms = {@(x) struct('L',{x}),@(x) struct('R',{x}),@(x) struct('in',{x}),@(x) struct('out',{x}),...
         @(x) struct('in',{x},'out',{x}),@(x) struct('L',{x},'R',{x}),@(x) struct('R',{x},'out',{x})};
for k=1:numel(forms)
    f = forms{k};
    out = cell(1,4);
    [out{:}] = piesos_quadopvar(prog1,f(v{1}),f(v{2}),f(v{3}),f(v{4}),v{5},v{6});
    assert_same_output(ref,out,sprintf('(b) variant %d',k));
end
assert(isequal(ref{3},ref{4}),'The same inputs must give ZopL=ZopR.');
% (c) which side is which: vector-valued state on the left only
[~,P,ZL,ZR] = piesos_quadopvar(prog1,struct('L',2,'R',1),{'s1'},[a b],1,1);
assert(isequal(Zdims(ZL),2*5) && isequal(Zdims(ZR),5) && isequal(size(P{1}),[10,5]),'dims.L and dims.R are mixed up.');
[~,P,ZL,ZR] = piesos_quadopvar(prog1,struct('out',2,'in',1),{'s1'},[a b],1,1);
assert(isequal(Zdims(ZL),2*5) && isequal(Zdims(ZR),5) && isequal(size(P{1}),[10,5]),'dims.out and dims.in are mixed up.');
end

function r = rows_LR(prog0,degs,opts)
% Row counts [rows of ZopL, rows of ZopR] for one scalar state in s1
[~,~,ZL,ZR] = piesos_quadopvar(prog0,1,{'s1'},[0.5 2],1,degs,opts);
r = [size(ZL,1),size(ZR,1)];
end

function case_degs_opts_forms()
% The forms of degs (array, cell, struct with in/out, struct with L and/or R)
% and of opts (top-level exclude/sep/otype, opts.L and/or opts.R)
rng(3);
a = 0.5;    b = 2;
prog0 = new_prog({'s1'});
call = @(degs,opts) piesos_quadopvar(prog0,1,{'s1'},[a b],1,degs,opts);
% degs: the same for L and R if an array, cell, struct in/out, or one of L, R
D = {2,struct('in',2,'out',2),struct('L',2),struct('R',2),struct('L',struct('in',2,'out',2)),...
     struct('R',struct('in',2,'out',2)),struct('L',2,'R',2),{{(0:2)'}},struct('L',{{{(0:2)'}}},'R',{{{(0:2)'}}})};
ref = cell(1,4);
[ref{:}] = call(D{1},struct());
for k=2:numel(D)
    out = cell(1,4);
    [out{:}] = call(D{k},struct());
    assert_same_output(ref,out,sprintf('degs form %d',k));
end
% degs: different for L and R
dref = struct('L',struct('in',1,'out',2),'R',struct('in',3,'out',0));
[ref{:}] = call(dref,struct());
assert(isequal(Zdims(ref{3}),nrows_3pi(1,false(1,3))) && isequal(Zdims(ref{4}),nrows_3pi(3,false(1,3))),...
    'degs.L.in and degs.R.in are mixed up.');
assert(isequal(size(ref{2}{1}),[5,9]),'Pmat has wrong size.');
out = cell(1,4);
[out{:}] = call(struct('L',struct('in',1,'out',2),'R',struct('in',{{3}},'out',{{0}})),struct());
assert_same_output(ref,out,'degs per-state cells');
% opts: exclude the multiplier (4 rows), the lower integral (3), or none (5)
exc = @(v) struct('exclude',logical(v(:)));
assert(isequal(rows_LR(prog0,1,struct()),[5,5]),'Default rows.');
assert(isequal(rows_LR(prog0,1,exc([1 0 0])),[4,4]),'Top-level exclude applies to both.');
assert(isequal(rows_LR(prog0,1,struct('L',exc([1 0 0]))),[4,4]),'Only opts.L: used for both.');
assert(isequal(rows_LR(prog0,1,struct('R',exc([0 1 0]))),[3,3]),'Only opts.R: used for both.');
assert(isequal(rows_LR(prog0,1,struct('L',exc([1 0 0]),'R',exc([0 1 0]))),[4,3]),'opts.L and opts.R.');
assert(isequal(rows_LR(prog0,1,struct('L',struct(),'R',exc([0 1 0]))),[5,3]),'An empty opts.L gives the defaults.');
assert(isequal(rows_LR(prog0,1,struct('R',exc([0 1 0]),'L',exc([0 0 1]),'pos',false,'sym',false)),[3,3]),...
    'opts.L and opts.R with flags.');
% opts.sep, on one side
dg = struct('in',2,'out',0);
assert(isequal(rows_LR(prog0,dg,struct('L',struct('sep',true))),[4,4]),'opts.L.sep applies to both.');
assert(isequal(rows_LR(prog0,dg,struct('L',struct(),'R',struct('sep',true))),[7,4]),'opts.R.sep only.');
end

function case_structure_1D()
% One variable, one state: different dummy degrees (rows of ZopL, ZopR) and
% primary degrees (monomials of ZsL, ZsR), so that a mix-up of left and right is detected
rng(4);
f3 = false(1,3);
for rr = {[1 2],[2 0],[0 3]}
    r = rr{1};
    for dd = {[1 2],[3 0]}
        d = dd{1};
        cc = cfg('rL',r(1),'rR',r(2),'oL',d(1),'oR',d(2));
        [prog,Pmat,ZL,ZR,prog0,spec] = run_structure(cc,'gen');
        check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
        mL = nrows_3pi(r(1),f3);        mR = nrows_3pi(r(2),f3);
        assert(isequal(size(Pmat{1}),[mL,mR]),'Pmat{1,1} must be %d x %d.',mL,mR);
        nd = mL*(d(1)+1)*mR*(d(2)+1);
        assert(prog.var.idx{2}-prog.var.idx{1}==nd,'Wrong number of decision variables.');
    end
end
end

function case_structure_blocks()
% Several blocks, a different number on the left and right, vector-valued
% states and own degrees per block: Pmat{i,j} couples the ith block row of
% ZopL and the jth block row of ZopR
rng(5);
a = 0.5;    b = 2;    f3 = false(1,3);
prog0 = new_prog({'s1'});
dims = struct('L',[1;2;1],'R',[2;1]);
dep = struct('L',[1;1;1],'R',[1;1]);
degs = struct('L',struct('in',{{1,2,0}},'out',{{1,0,2}}),'R',struct('in',{{0,1}},'out',{{1,1}}));
[prog,Pmat,ZL,ZR] = piesos_quadopvar(prog0,dims,{'s1'},[a b],dep,degs,struct());
mLx = [nrows_3pi(1,f3); 2*nrows_3pi(2,f3); nrows_3pi(0,f3)];
mRx = [2*nrows_3pi(0,f3); nrows_3pi(1,f3)];
check_Zop(ZL,mLx,{{(0:1)'},{(0:2)'},{0}},{'s1'},{'s1'});
check_Zop(ZR,mRx,{{0},{(0:1)'}},{'s1'},{'s1'});
spec = struct('vused',{{'s1'}},'dom',[a b],'EL',{{{(0:1)'},{0},{(0:2)'}}},'ER',{{{(0:1)'},{(0:1)'}}},...
    'g',{no_multiplier()},'mode','gen');
check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
end

function case_structure_spaces()
% Left and right states in different variables: both map to the space of all
% variables on which some state depends (left variables first). Mixed
% (nested), reversed, disjoint, three variables, an unused variable, and
% variables in different order
rng(6);
a = 0.5;    b = 2;    dm2 = [-1 0.7];     dm3 = [1 3];
cs = {};
% nested: left x(s1), right y(s1,s2)
cs{end+1} = cfg('vL',{'s1'},'vR',{'s1','s2'},'mL',[a b],'mR',[a b;dm2],'pL',1,'pR',[1 1],...
                'rL',1,'rR',[1 2],'oL',[1 1],'oR',[2 1]);
% reversed: left y(s1,s2), right x(s2)
cs{end+1} = cfg('vL',{'s1','s2'},'vR',{'s2'},'mL',[a b;dm2],'mR',dm2,'pL',[1 1],'pR',1,...
                'rL',[2 1],'rR',1,'oL',[1 1],'oR',[1 1]);
% disjoint: left x(s1), right y(s2)
cs{end+1} = cfg('vL',{'s1'},'vR',{'s2'},'mL',[a b],'mR',dm2,'rL',1,'rR',2,'oL',[1 0],'oR',[0 2]);
% disjoint, reversed order of the primary variables: (s2,s1)
cs{end+1} = cfg('vL',{'s2'},'vR',{'s1'},'mL',dm2,'mR',[a b],'rL',1,'rR',2,'oL',[1 0],'oR',[0 2]);
% three variables, a shared one: (s1,s3) on the left, (s2,s3) on the right
cs{end+1} = cfg('vL',{'s1','s3'},'vR',{'s2','s3'},'mL',[a b;dm3],'mR',[dm2;dm3],'pL',[1 1],'pR',[1 1],...
                'rL',[1 0],'rR',[1 1],'oL',[1 0 0],'oR',[0 1 1]);
% an unused variable in the left list: s4 does not appear
cs{end+1} = cfg('vL',{'s1','s4'},'vR',{'s1'},'mL',[a b;0 1],'mR',[a b],'pL',[1 0],'pR',1,...
                'rL',[1 0],'rR',2,'oL',[1 0],'oR',[2 0]);
% blocks of different dependence on both sides, one state without variables
cs{end+1} = cfg('dL',[1;1;1],'dR',2,'vL',{'s1','s2'},'vR',{'s2'},'mL',[a b;dm2],'mR',dm2,...
                'pL',[0 0;1 0;1 1],'pR',1,'rL',[1 1],'rR',1,'oL',[1 0],'oR',[1 1]);
% variables in different order on the left and right (dummy degree 0)
cs{end+1} = cfg('vL',{'s2','s1'},'vR',{'s1','s2'},'mL',[dm2;a b],'mR',[a b;dm2],'pL',[1 1],'pR',[1 1],...
                'rL',[0 0],'rR',[0 0],'oL',[2 1],'oR',[1 2]);
for c=1:numel(cs)
    cc = cs{c};
    [prog,Pmat,ZL,ZR,prog0,spec] = run_structure(cc,'gen');
    check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
    for k=1:numel(Pmat)
        assert(all(ismember(Pmat{k}.varname,spec.vused)),'Case %d: Pmat depends on an unused variable.',c);
    end
end
end

function case_swap()
% Exchanging the left and right inputs exchanges ZopL and ZopR, and
% transposes the block structure of Pmat
rng(7);
a = 0.5;    b = 2;
dims = struct('L',[1;2],'R',1);
vars = struct('L',{{'s1'}},'R',{{'s1','s2'}});
dom = struct('L',[a b],'R',[a b;-1 0.7]);
dep = struct('L',[1;1],'R',[1 1]);
degs = struct('L',struct('in',{{1,2}},'out',{{1,0}}),'R',struct('in',[1 2],'out',[0 1]));
sw = @(s) struct('L',{s.R},'R',{s.L});
prog0 = new_prog({'s1','s2'});
[p1,P1,ZL1,ZR1] = piesos_quadopvar(prog0,dims,vars,dom,dep,degs,struct());
[p2,P2,ZL2,ZR2] = piesos_quadopvar(prog0,sw(dims),sw(vars),sw(dom),sw(dep),sw(degs),struct());
assert(isequal(ZL1,ZR2) && isequal(ZR1,ZL2),'Exchanging L and R must exchange ZopL and ZopR.');
assert(isequal(size(P1),fliplr(size(P2))),'Exchanging L and R must transpose the block structure.');
for i=1:size(P1,1)
    assert(isequal(size(P1{i,1}),fliplr(size(P2{1,i}))),'Block (%d,1) must be transposed.',i);
end
assert(p1.var.idx{2}-p1.var.idx{1}==p2.var.idx{2}-p2.var.idx{1},'Number of decision variables must agree.');
end

function case_modes_same_dims()
% ZopL~=ZopR with equal output dimensions and primary degrees: Pmat general,
% symmetric or PSD
rng(8);
a = 0.5;    b = 2;    dm2 = [-1 0.7];
cs = {};
% exclude the lower integral on the left, the upper on the right (3 rows each)
cs{end+1} = cfg('rL',1,'rR',1,'oL',1,'oR',1,'xL',[0 1 0],'xR',[0 0 1]);
% the same on two blocks, one vector-valued
cs{end+1} = cfg('dL',[1;2],'dR',[1;2],'pL',[1;1],'pR',[1;1],'rL',1,'rR',1,'oL',1,'oR',1,'xL',[0 1 0],'xR',[0 0 1]);
% disjoint variables, 5 rows each
cs{end+1} = cfg('vL',{'s1'},'vR',{'s2'},'mL',[a b],'mR',dm2,'rL',1,'rR',1,'oL',[1 1],'oR',[1 1]);
% both in (s1,s2): different dummy degrees, 15 rows each
cs{end+1} = cfg('vL',{'s1','s2'},'vR',{'s1','s2'},'mL',[a b;dm2],'mR',[a b;dm2],'pL',[1 1],'pR',[1 1],...
                'rL',[1 0],'rR',[0 1],'oL',[1 0],'oR',[1 0]);
% nested spaces: x(s1) on the left (dummy degree 3), y(s1,s2) on the right
% (degree 0 in both): 9 rows each
cs{end+1} = cfg('vL',{'s1'},'vR',{'s1','s2'},'mL',[a b],'mR',[a b;dm2],'pL',1,'pR',[1 1],...
                'rL',3,'rR',[0 0],'oL',[1 0],'oR',[1 0]);
for c=1:numel(cs)
    for md = {'gen','sym','pos'}
        [prog,Pmat,ZL,ZR,prog0,spec] = run_structure(cs{c},md{1});
        check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
        assert(~isequal(ZL,ZR),'Case %d: ZopL and ZopR must differ.',c);
    end
end
end

function case_psatz_LR()
% psatz multipliers with ZopL~=ZopR (opts.pos): all terms, orders and logical
% selections; monomials cut by ceil(qk/2)
rng(9);
a = 0.5;    b = 2;    dm2 = [-1 0.7];
% 1D, exclude lower on the left and upper on the right, orders 1,2
for q=1:2
    cc = cfg('rL',1,'rR',1,'oL',4,'oR',4,'xL',[0 1 0],'xR',[0 0 1]);
    [prog,Pmat,ZL,ZR,prog0,spec] = run_structure(cc,'pos',struct('psatz',q));
    spec.g = psatz_terms([a b],q,1:q+1);
    spec.EL = cut_E(spec.EL,q);         spec.ER = cut_E(spec.ER,q);
    check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
    assert(prog.var.num==q+1,'Expected %d SOS variables.',q+1);
end
% logical selection of the terms, 1D
for sel = {[true;false;true],[true false true]}
    cc = cfg('rL',0,'rR',0,'oL',3,'oR',3,'xL',[0 1 0],'xR',[0 0 1]);
    [prog,Pmat,ZL,ZR,prog0,spec] = run_structure(cc,'pos',struct('psatz',sel{1}));
    spec.g = {@(s) (b-s(1))^2,@(s) (s(1)-a)^2};
    spec.EL = cut_E(spec.EL,2);         spec.ER = cut_E(spec.ER,2);
    check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
end
% disjoint variables (s1) on the left, (s2) on the right: psatz in both, orders [1 2]
cc = cfg('vL',{'s1'},'vR',{'s2'},'mL',[a b],'mR',dm2,'rL',1,'rR',1,'oL',[2 3],'oR',[2 3]);
[prog,Pmat,ZL,ZR,prog0,spec,~,dout] = run_structure(cc,'pos',struct('psatz',[1 2]));
spec.g = psatz_terms(dout,[1 2],1:6);
spec.EL = cut_E(spec.EL,[1 2]);         spec.ER = cut_E(spec.ER,[1 2]);
check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
assert(prog.var.num==6,'Expected 6 SOS variables.');
% nested (x(s1) left, y(s1,s2) right, 9 rows each), scalar psatz: the same
% order in the variables used
cc = cfg('vL',{'s1'},'vR',{'s1','s2'},'mL',[a b],'mR',[a b;dm2],'pL',1,'pR',[1 1],...
         'rL',3,'rR',[0 0],'oL',[2 2],'oR',[2 2]);
[prog,Pmat,ZL,ZR,prog0,spec,~,dout] = run_structure(cc,'pos',struct('psatz',1));
spec.g = psatz_terms(dout,[1 1],1:4);
spec.EL = cut_E(spec.EL,[1 1]);         spec.ER = cut_E(spec.ER,[1 1]);
check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
% an unused primary variable: a scalar psatz is applied to the used one only
cc = cfg('vL',{'s1','s2'},'vR',{'s1'},'mL',[a b;dm2],'mR',[a b],'pL',[1 0],'pR',1,...
         'rL',[1 0],'rR',1,'oL',[2 0],'oR',[2 0]);
[prog,Pmat,ZL,ZR,prog0,spec] = run_structure(cc,'pos',struct('psatz',1));
spec.g = psatz_terms([a b],1,1:2);
spec.EL = cut_E(spec.EL,1);             spec.ER = cut_E(spec.ER,1);
check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
assert(~ismember('s2',Pmat{1}.varname),'Pmat must not depend on the unused variable s2.');
end

function case_sequential()
% Two declarations in the same program, one general and one positive with a
% psatz multiplier: the decision variables must not overlap
rng(12);
a = 0.5;    b = 2;
prog0 = new_prog({'s1'});
degs1 = struct('L',struct('in',1,'out',1),'R',struct('in',2,'out',2));
[prog1,P1,ZL1,ZR1] = piesos_quadopvar(prog0,1,{'s1'},[a b],1,degs1,struct());
[prog2,P2,ZL2,ZR2] = piesos_quadopvar(prog1,1,{'s1'},[a b],1,struct('in',1,'out',3),struct('pos',true,'psatz',1));
assert(prog1.var.num==1 && prog2.var.num==3,'Expected 1+2 declared variables.');
spec1 = struct('vused',{{'s1'}},'dom',[a b],'EL',{{{(0:1)'}}},'ER',{{{(0:2)'}}},'g',{no_multiplier()},'mode','gen');
spec2 = struct('vused',{{'s1'}},'dom',[a b],'EL',{{{(0:2)'}}},'ER',{{{(0:2)'}}},'g',{{@(s) b-s(1),@(s) s(1)-a}},'mode','pos');
check_sos_LR(prog1,P1,ZL1,ZR1,spec1,prog0);
check_sos_LR(prog2,P2,ZL2,ZR2,spec2,prog1);
assert(isempty(intersect(P1{1}.dvarname,P2{1}.dvarname)),'The two matrices must not share decision variables.');
end

function case_no_variables()
% No spatial variables: states in R^2 and R^1, ZopL, ZopR the identity, Pmat
% a constant matrix
rng(13);
prog0 = new_prog({'s1'});
[prog,Pmat,ZL,ZR] = piesos_quadopvar(prog0,struct('L',[2;1],'R',1),{},zeros(0,2),struct('L',zeros(2,0),'R',zeros(1,0)));
assert(isequal(Zdims(ZL),[2;1]) && isequal(Zdims(ZR),1),'ZopL and ZopR must be the identity on the states.');
spec = struct('vused',{cell(1,0)},'dom',zeros(0,2),'EL',{{cell(1,0),cell(1,0)}},'ER',{{cell(1,0)}},'g',{no_multiplier()},'mode','gen');
check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
[prog,Pmat,ZL,ZR] = piesos_quadopvar(prog0,[2;1],{},zeros(0,2),zeros(2,0),0,struct('pos',true));
spec = struct('vused',{cell(1,0)},'dom',zeros(0,2),'EL',{{cell(1,0),cell(1,0)}},'ER',{{cell(1,0),cell(1,0)}},'g',{no_multiplier()},'mode','pos');
check_sos_LR(prog,Pmat,ZL,ZR,spec,prog0);
end

function case_errors()
prog = new_prog({'s1'});
dg = struct('in',1,'out',1);
q1 = @(dims,degs,opts) piesos_quadopvar(prog,dims,{'s1'},[0 1],1,degs,opts);
q2 = @(dims,dep,degs,opts) piesos_quadopvar(prog,dims,{'s1'},[0 1],dep,degs,opts);
degsLR = struct('L',struct('in',1,'out',1),'R',struct('in',2,'out',1));
% symmetric and positive: block dimensions and primary monomials must match
assert_err(@() q1(1,degsLR,struct('pos',true)),'Output dimensions');
assert_err(@() q1(1,degsLR,struct('sym',true)),'Output dimensions');
assert_err(@() q2(struct('L',[1;1],'R',1),struct('L',[1;1],'R',1),dg,struct('pos',true)),'Output dimensions');
assert_err(@() q2(struct('L',1,'R',[1;1]),struct('L',1,'R',[1;1]),dg,struct('sym',true)),'Output dimensions');
assert_err(@() q1(struct('L',1,'R',2),dg,struct('sym',true)),'Output dimensions');
degsO = struct('L',struct('in',1,'out',1),'R',struct('in',1,'out',2));
assert_err(@() q1(1,degsO,struct('pos',true)),'Degrees of the primary');
assert_err(@() q1(1,degsO,struct('sym',true)),'Degrees of the primary');
q1(1,degsO,struct());                                   % allowed if not symmetric
q1(1,degsLR,struct());
q2(struct('L',[1;1],'R',1),struct('L',[1;1],'R',1),dg,struct());
q1(1,dg,struct('psatz',0));                             % order 0: no multiplier
% the same name must have the same domain
assert_err(@() piesos_quadopvar(prog,1,{'s1'},struct('L',[0 1],'R',[0 2]),1,dg,struct()),'same domain');
% the inputs for left and right
assert_err(@() piesos_quadopvar(prog,struct('x',1),{'s1'},[0 1],1,dg,struct()),'array, or as a struct');
assert_err(@() piesos_quadopvar(prog,struct('L',1,'out',1),{'s1'},[0 1],1,dg,struct()),'only once');
assert_err(@() piesos_quadopvar(prog,1,struct('R',{{'s1'}},'in',{{'s1'}}),[0 1],1,dg,struct()),'only once');
assert_err(@() piesos_quadopvar(prog,1,5,[0 1],1,dg,struct()),'cellstr');
assert_err(@() piesos_quadopvar(prog,1,{'s1','s1'},[0 1;0 1],[1 1],dg,struct()),'distinct');
assert_err(@() piesos_quadopvar(prog,1,{'s1'},[0 1;0 1],1,dg,struct()),'Domains');
assert_err(@() piesos_quadopvar(prog,struct('L',1,'R',[1;1]),{'s1'},[0 1],1,dg,struct()),'Dependency array');
assert_err(@() piesos_quadopvar(prog,1,{'s1'},[0 1],[1 1],dg,struct()),'Dependency array');
assert_err(@() piesos_quadopvar(prog,1,{'s1'},[0 1],struct('L',1,'R',[1 1]),dg,struct()),'Dependency array');
% options
assert_err(@() q1(1,dg,5),'struct');
assert_err(@() q1(1,dg,struct('L',5)),'must be specified as');
assert_err(@() q1(1,dg,struct('R',5,'L',struct())),'must be specified as');
assert_err(@() q1(1,dg,struct('otype',[2 2])),'currently not supported');
assert_err(@() q1(1,dg,struct('L',struct('otype',[2 2]))),'currently not supported');
assert_err(@() q1(1,dg,struct('L',struct(),'R',struct('otype',[2 2]))),'currently not supported');
% invalid degrees and psatz, as for piesos_poslpivar
assert_err(@() q1(1,struct('in',1,'out',-1),struct()),'nonnegative integers');
assert_err(@() q1(1,struct('in',1,'out',1.5),struct()),'nonnegative integers');
assert_err(@() q1(1,struct('in',1),struct()),'fields ''in'' and ''out''');
assert_err(@() q1(1,struct('L',struct('in',1)),struct()),'fields ''in'' and ''out''');
assert_err(@() q1(1,struct('in',1,'out',[1 2]),struct()),'1 x M');
assert_err(@() q1(1,dg,struct('pos',true,'psatz',-1)),'nonnegative integers');
assert_err(@() q1(1,dg,struct('pos',true,'psatz',[1 2])),'1 x N array');
assert_err(@() q1(1,dg,struct('pos',true,'psatz','a')),'1 x N array');
assert_err(@() q1(1,dg,struct('pos',true,'psatz',[false;false])),'must be true');
prog2 = new_prog({'s1','s2'});
assert_err(@() piesos_quadopvar(prog2,1,{'s1','s2'},[0 1;0 1],[1 0],1,struct('pos',true,'psatz',[1 1])),'no state depends');
% a variable that only the right operator uses is used for the psatz too
piesos_quadopvar(prog2,1,struct('L',{{'s1'}},'R',{{'s2'}}),[0 1],1,1,struct('pos',true,'psatz',[1 1]));
end
