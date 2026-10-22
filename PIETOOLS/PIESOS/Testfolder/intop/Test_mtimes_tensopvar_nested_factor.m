% Test_mtimes_tensopvar_nested_factor.m -- Tests the composition of an
% 'intop' with a 'tensopvar' whose factors are ITSELF nested 'tensopvar'
% objects (see tensopvar.m and @tensopvar/initialize.m), as in the
% nonlinear terms of a PIE: the pointwise product (otype=[1,2]) of
% operators applied to different copies of the state is a factor in a
% Kronecker product of operators (otype=[2,2]). The composition is a
% functional on the tensor product of the inputs of all the leaf operators.
%
% Ground truth, independent of mtimes_tensopvar: the functional Aop applied
% to the output of Bop, evaluated by recursion over the nested operators
% (apply_nested below: Kronecker or pointwise product of the outputs of the
% factors, the dummy variable of an output variable s numbered by the
% factors that precede it, s_1, s_2, ... as in apply_op), vs. the
% composition Cop applied to the tensor product of the polynomial states of
% the leaf operators, each in its own dummy variables, as in
% Test_intop_mtimes_tensopvar.m.
%
% Cases: 1 and 2 spatial variables, the pointwise product of 2 and 3
% operators, the sibling operator as a 'sopvar' or a (single-factor)
% 'tensopvar', before or after the nested operator, terms in the nested
% operator and in the outer operator, 'dpvar' kernels (with different
% decision variables in each region), pointwise product of products, nested
% operators in nested operators, and variables that only some of the
% operators map from.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/21/2026: Names of the variable groups S1 (input only) and S2
%   (output only) swapped, to match 'sopvar' and mtimes_tensopvar;
%
% NOTES
% - Nestings that mtimes_tensopvar does not support (Kronecker product as a
%   factor in a pointwise product) are only checked to be rejected with a
%   'not supported' message, section (iv).
% - Sparse random operators can give a zero reference value: check_deep
%   draws again (up to 20 times, error id
%   'Test_mtimes_tensopvar_nested_factor:vacuous'); a zero reference is
%   never accepted.
% - All kernels of an 'intop' share one list Zd: thin_dvars gives each
%   kernel a different subset of Zd by zeroing the rows of the others.

tol = 1e-8;
ncases = 0;

%% (i) pointwise product of operators as a factor in a Kronecker product,
%   B(x1,x2,x3) = N(x1,x2) o S(x3),   N(x1,x2)(s) = (P1 x1)(s)*(P2 x2)(s) [*(P3 x3)(s)]
%   Ntot spatial variables, all operators mapping from and to all of them, unless 'dep'
%   specifies the variables (rows: P1,...,P3 and S) that they map from and to. Fields of each
%   case: nleaves (operators in N), order of N and S in B, sibt (S is a tensopvar, not a
%   sopvar), nNterms (terms in N), nBterms (terms in B), type (otype of B), q (decision
%   variables in the kernels of Aop), dnsty (density of the operators), zero_region
cases = {
    struct('Ntot',1,'nleaves',2,'order','NS','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',0)
    struct('Ntot',1,'nleaves',2,'order','SN','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',0)
    struct('Ntot',1,'nleaves',2,'order','NS','sibt',true, 'nNterms',2,'nBterms',1,'type',[2,2],'q',0)
    struct('Ntot',1,'nleaves',2,'order','SN','sibt',true, 'nNterms',1,'nBterms',1,'type',[2,2],'q',0)
    struct('Ntot',1,'nleaves',3,'order','NS','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',0,'dnsty',0.3)
    struct('Ntot',1,'nleaves',3,'order','SN','sibt',true, 'nNterms',1,'nBterms',1,'type',[2,2],'q',0,'dnsty',0.3)
    struct('Ntot',1,'nleaves',2,'order','NS','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',5)
    struct('Ntot',1,'nleaves',2,'order','NS','sibt',true, 'nNterms',2,'nBterms',2,'type',[2,2],'q',0)
    struct('Ntot',1,'nleaves',2,'order','SN','sibt',true, 'nNterms',1,'nBterms',2,'type',[2,2],'q',0)
    struct('Ntot',1,'nleaves',2,'order','NS','sibt',true, 'nNterms',1,'nBterms',1,'type',[1,2],'q',0)
    % an integral that is zero in the nested operator of the second term, but not of the first
    struct('Ntot',1,'nleaves',2,'order','NS','sibt',true, 'nNterms',1,'nBterms',2,'type',[2,2],'q',0,'zero_region',true)
    struct('Ntot',1,'nleaves',2,'order','SN','sibt',false,'nNterms',2,'nBterms',2,'type',[2,2],'q',0,'zero_region',true)
    % two spatial variables: operators map from and to different combinations of them
    struct('Ntot',2,'nleaves',2,'order','NS','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',0,...
           'dep',struct('in',[1 1;1 0;0 1],'out',[1 1;1 0;0 1]))
    struct('Ntot',2,'nleaves',2,'order','SN','sibt',true, 'nNterms',1,'nBterms',1,'type',[2,2],'q',0,...
           'dep',struct('in',[1 0;1 1;1 1],'out',[1 0;1 1;0 1]))
    struct('Ntot',2,'nleaves',2,'order','NS','sibt',true, 'nNterms',1,'nBterms',1,'type',[2,2],'q',0,...
           'dep',struct('in',[1 0;1 0;0 1],'out',[1 0;1 0;0 1]))
    % a variable that only the first operator maps from (integrated out): S1 variable % DJ (CC), 10/21/2026
    struct('Ntot',2,'nleaves',2,'order','NS','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',0,...
           'dep',struct('in',[1 1;1 0;1 0],'out',[1 0;1 0;1 0]))
    struct('Ntot',2,'nleaves',2,'order','SN','sibt',false,'nNterms',1,'nBterms',1,'type',[2,2],'q',0,...
           'dep',struct('in',[1 0;1 1;1 0],'out',[1 0;1 0;1 0]))
    };
for c=1:numel(cases)
    rng(1000+c);
    check_case(cases{c},tol,sprintf('(i) case %d',c));
    ncases = ncases+1;
end

%% (ii) nested operators in nested operators: B = {{N1 o S1} o N2} or {N1 o {S1 o N2}}
for c = 1:2
    rng(2000+c);
    check_deep(c,tol,sprintf('(ii) case %d',c));
    ncases = ncases+1;
end

%% (iii) nested Kronecker product as the only factor, without variables that are passed through
dims_in = [1,1; 1,1];
vars_in = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
dom_in = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
depmat_in = struct('in',logical(eye(2)),'out',logical(eye(2)));
rng(7301);
Rinner = rand_tensopvar(dims_in,vars_in,dom_in,depmat_in,1,[2,2],1,0.5);
Bop = tensopvar();      Bop.otype = [2,2];      Bop.ops = {Rinner};     Bop = initialize(Bop);
pdegs = ones(1,numel(Bop.vars.out));
Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);
pv = polynomial(Bop.vars.in(:));
leaves = {rand_poly([1,1],pv(1),2,10),rand_poly([1,1],pv(2),2,10)};
ABv = apply_op(Aop,apply_nested(Bop,leaves,1,zeros(1,numel(Bop.vars.out)),Bop.vars.out));
Cop = mtimes_tensopvar(Aop,Bop);
Cv = apply_op(Cop,kron(leaves{1},leaves{2}));
compare(ABv,Cv,tol,'(iii) nested Kronecker product');
ncases = ncases+1;

%% (iv) unsupported nestings must still be rejected, with an informative message
% Kronecker product of operators as a factor in a pointwise product
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
dep2 = struct('in',true(2,1),'out',true(2,1));
dep1 = struct('in',true,'out',true);
rng(5000);
Nk = rand_tensopvar([1,1;1,1],vars,dom,dep2,1,[2,2],1,0.5);
Sk = rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],1,0.5);
Bk = tensopvar();       Bk.otype = [1,2];       Bk.ops = {Nk,Sk};       Bk = initialize(Bk);
Ak = rand_intop([1,Bk.matdim(1)],Bk.vars.out,Bk.dom.out,Bk.deps.out,1);
expect_error(@() mtimes_tensopvar(Ak,Bk),'not supported');

fprintf('Test_mtimes_tensopvar_nested_factor passed (%d cases).\n',ncases)


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function check_case(cfg,tol,label)
% Build the nested operator B and the functional A, and compare the composition with the ground truth
Ntot = cfg.Ntot;
svars = arrayfun(@(i) ['s',num2str(i)],1:Ntot,'UniformOutput',false);
nl = cfg.nleaves;
if isfield(cfg,'dep')
    depin = logical(cfg.dep.in);        depout = logical(cfg.dep.out);
else
    depin = true(nl+1,Ntot);            depout = true(nl+1,Ntot);
end
% Variables of the operators: those that some operator maps from or to
is_in = any(depin,1);                   is_out = any(depout,1);
vars = struct('in',{svars(is_in)},'out',{svars(is_out)});
depin = depin(:,is_in);                 depout = depout(:,is_out);
N = numel(vars.in);                     M = numel(vars.out);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],M,1));
m = randi(2);                   % row dimension of the pointwise product
nin = randi(2,nl,1);            % column dimensions of the operators in the pointwise product
n3 = randi(2);                  % column dimension of the sibling
depN = struct('in',depin(1:nl,:),'out',depout(1:nl,:));
depS = struct('in',depin(nl+1,:),'out',depout(nl+1,:));
pdegs = struct('in',ones(1,N),'out',ones(1,M));
dnsty = 0.5;
if isfield(cfg,'dnsty')
    dnsty = cfg.dnsty;
end
if cfg.type(1)==1
    dims_S = [m,n3];            % elementwise product with the sibling
else
    dims_S = [randi(2),n3];
end
Sop = rand_tensopvar(dims_S,vars,dom,depS,1,[2,2],pdegs,dnsty);
if ~cfg.sibt
    Sop = Sop.ops{1,1};         % single 'sopvar' factor
end
Nterm = cell(cfg.nBterms,1);
for it=1:cfg.nBterms
    Nterm{it} = rand_tensopvar([m*ones(nl,1),nin],vars,dom,depN,cfg.nNterms,[1,2],pdegs,dnsty);
end
if isfield(cfg,'zero_region') && cfg.zero_region
    % Zero the upper integral in the first operator of the nested operator in the second term
    op = Nterm{2}.ops{1,1};
    P = op.params;
    P{end} = sparse(size(P{end},1),size(P{end},2));
    op.params = P;
    Nterm{2}.ops{1,1} = op;
end
Bop = tensopvar();      Bop.otype = cfg.type;
if strcmp(cfg.order,'NS')
    Bop.ops = [Nterm,repmat({Sop},cfg.nBterms,1)];
else
    Bop.ops = [repmat({Sop},cfg.nBterms,1),Nterm];
end
Bop = initialize(Bop);

% States of the leaf operators, in the order in which they appear in B
pvars = polynomial(Bop.vars.in(:));
xl = cell(1,nl);
for l=1:nl
    xl{l} = rand_poly([nin(l),1],pvars(depN.in(l,:)),2,10);
end
x3 = rand_poly([n3,1],pvars(depS.in(1,:)),2,10);
if strcmp(cfg.order,'NS')
    leaves = [xl,{x3}];     deps_in = [depN.in;depS.in];
else
    leaves = [{x3},xl];     deps_in = [depS.in;depN.in];
end
compose_and_compare(Bop,leaves,deps_in,cfg.q,tol,label);
end

function check_deep(variant,tol,label)
% Nested operators in nested operators, in one spatial variable, with sparse operators: the
% reference value may be zero, in which case the operators are drawn again
for attempt = 1:20
    try
        check_deep_once(variant,tol,label);
        return
    catch ME
        if ~strcmp(ME.identifier,'Test_mtimes_tensopvar_nested_factor:vacuous')
            rethrow(ME)
        end
    end
end
error("Test_mtimes_tensopvar_nested_factor (%s): the reference value is zero in each of the draws.",label)
end

function check_deep_once(variant,tol,label)
% Operators: multipliers, except for one integral operator (S1 in variant 1, N2(2) in variant 2),
% to keep the number of integrals low
vars = struct('in',{{'s'}},'out',{{'s'}});
dom = struct('in',[0,1],'out',[0,1]);
dep1 = struct('in',true,'out',true);
pd = struct('in',1,'out',1);
pvar s
mult = @() mat2tensopvar(randn+randn*s+randn*s^2,vars,dom).ops{1,1};      % multiplier operator
intl = @() rand_tensopvar([1,1],vars,dom,dep1,1,[2,2],pd,0.5).ops{1,1};    % random 3-PI operator
pw = @(A,B) pointwise(A,B);
if variant==1
    % {{N1 o S1} o N2}, N1 = P1*P2, N2 = P3*P4
    N1 = pw(mult(),mult());     S1 = intl();      N2 = pw(mult(),mult());
    R12 = tensopvar();      R12.otype = [2,2];      R12.ops = {N1,S1};      R12 = initialize(R12);
    Bop = tensopvar();      Bop.otype = [2,2];      Bop.ops = {R12,N2};
else
    % {N1 o {S1 o N2}}, N1 = P1*P2, N2 = P3*P4, with one integral operator in N2
    N1 = pw(mult(),mult());     S1 = mult();      N2 = pw(mult(),intl());
    R12 = tensopvar();      R12.otype = [2,2];      R12.ops = {S1,N2};      R12 = initialize(R12);
    Bop = tensopvar();      Bop.otype = [2,2];      Bop.ops = {N1,R12};
end
Bop = initialize(Bop);
nleaves = 5;
pvars = polynomial(Bop.vars.in(:));
leaves = cell(1,nleaves);
for l=1:nleaves
    leaves{l} = rand_poly([1,1],pvars,2,10);
end
compose_and_compare(Bop,leaves,true(nleaves,1),0,tol,label);
end

function N = pointwise(A,B)
% Pointwise product of the output of the two operators
N = tensopvar();        N.otype = [1,2];        N.ops = {A,B};      N = initialize(N);
end
function compose_and_compare(Bop,leaves,deps_in,q,tol,label)
% Compare Aop*Bop with the functional Aop applied to Bop applied to the states
M = numel(Bop.vars.out);        N = numel(Bop.vars.in);
nr = randi(2);
Aop = rand_intop([nr,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,ones(1,M),q);
if q>0
    Aop = thin_dvars(Aop);
    d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:q,'UniformOutput',false),num2cell(randn(q,1)));
    Anum = eval_intop_dpvar(Aop,d0map);
else
    Anum = Aop;
end
% Ground truth: the functional applied to the output of Bop
ABv = apply_op(Anum,apply_nested(Bop,leaves,1,zeros(1,M),Bop.vars.out));
% Composition, acting on the tensor product of the states, each in its own dummy variables
Cop = mtimes_tensopvar(Aop,Bop);
pvars = polynomial(Bop.vars.in(:));
slot = zeros(1,N);      v_full = 1;
for l=1:numel(leaves)
    vl = leaves{l};
    for k=find(deps_in(l,:))
        slot(k) = slot(k)+1;
        vl = subs(vl,pvars(k),polynomial({[Bop.vars.in{k},'_dum',num2str(slot(k))]}));
    end
    v_full = kron(v_full,vl);
end
Cv = apply_op(Cop,v_full);
if q>0
    Cv = eval_dpvar_at(Cv,d0map);
end
compare(ABv,Cv,tol,label);
end

function [Y,pos,nslots] = apply_nested(R,leaves,pos,off,vars_out)
% Apply the (nested) tensor-PI operator R to the states in leaves(pos:end), one for each of
% the 'sopvar' operators in the first term of R, in order. The output variable vars_out{k}
% is replaced by vars_out{k}_(off(k)+1), vars_out{k}_(off(k)+2),... for the factors in a
% Kronecker product of the outputs (as in apply_intop), whereas the output of the factors in
% a pointwise product share the same variable.
%
% OUTPUTS
% - Y:      polynomial array, the output of R;
% - pos:    index of the next state in leaves;
% - nslots: 1 x M array, the number of variables s_1,s_2,... for each output variable;
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
    % Each term acts on the same states
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

function compare(ABv,Cv,tol,label)
err = cleanpoly(ABv-Cv,1e-10);
nrmref = full(max(max(abs(double(ABv)))));
if nrmref<=1e-6
    error('Test_mtimes_tensopvar_nested_factor:vacuous',"Test_mtimes_tensopvar_nested_factor (%s): the reference value is zero, the test is vacuous.",label)
elseif ~isdouble(err) || max(max(abs(double(err))))>tol*max(1,nrmref)
    error("Test_mtimes_tensopvar_nested_factor failed (%s): residual %.3e exceeds tolerance.",label,max(max(abs(double(err)))))
end
end

function Aop = thin_dvars(Aop)
% Retain a different random subset of the decision variables in each parameter of Aop
% (the rows of the others, in the blocks of rows of each row of the matrix, are set to zero)
nd = numel(Aop.Zd);
mm = size(Aop,1);
for k=1:numel(Aop.params)
    P = Aop.params{k};
    if nd==0,       continue,       end
    J = sort(randperm(nd,ceil(nd/2)));
    rows = reshape(((0:mm-1)'*(nd+1) + [1,1+J]).',[],1);
    P(setdiff(1:size(P,1),rows),:) = 0;
    Aop.params{k} = P;
end
end

function expect_error(fun,str)
try
    fun();
catch ME
    if ~contains(ME.message,str)
        error("Test_mtimes_tensopvar_nested_factor failed: unexpected error message '%s'.",ME.message)
    end
    return
end
error("Test_mtimes_tensopvar_nested_factor failed: expected an error, but none was thrown.")
end
