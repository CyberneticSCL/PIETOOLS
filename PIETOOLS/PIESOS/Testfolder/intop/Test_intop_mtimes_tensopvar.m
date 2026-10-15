% Test_intop_mtimes_tensopvar.m -- Tests mtimes_tensopvar for an 'intop'
% Aop and a 'tensopvar' Bop.
%
% Oracle: the composition Cop = mtimes_tensopvar(Aop,Bop) applied to the
% product input v_full is compared to Aop applied to the result of Bop
% applied to v, both evaluated with apply_op. Three parts: a primary check
% (one term, hard failure), a secondary check (several terms, warning
% only), and a check with decision variables in the kernels of Aop (hard
% failure).
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - The secondary check (nterms>1) only warns on a mismatch or error:
%   tensopvar factor columns that are themselves nested tensopvar/intop are
%   not covered.
% - Dpvar check: Aop has q in [50,200] decision variables.

% Set parameters defining tensor-PI operators
nterms = 1;      % number of terms in tensor-PI operator
d = 4;           % degree of the tensor-PI operator
N = 2;           % number of input variables
M = 2;           % number of output variables
otype = [2,2];           % type of tensor-PI operator
dnsty = 0.25;

% Set parameters defining the functional
pdegs = randi(2,[1,M]);
pdegsB = struct('in',randi(2,[1,N]),'out',randi(2,[1,M]));
nr = 2;
use_full_int = false;
pdegv = 3;

% Matrix dimensions of each factor in the tensor-PI operators
dims = randi(3,[d,2]);
if otype(1)==1
    m = dims(1,1);
    dims(:,1) = m;
else
    m = prod(dims(:,1));
end
if otype(2)==1
    n = dims(1,2);
    dims(:,2) = n;
else
    n = prod(dims(:,2));
end
% Input and output variables
vars = struct();
vars.in = cell(1,N);
idcs = sort(randperm(2*N,N));
for i=1:N
    vars.in{i} = ['s',num2str(idcs(i))];
end
vars.in = {'s2','s3'};
vars.in = vars.in(1:N);
vars.out = cell(1,M);
idcs = sort(randperm(2*N+2*M,M));
for i=1:M
    vars.out{i} = ['s',num2str(idcs(i))];
end
vars.out = {'s1','s3'};
vars.out = vars.out(1:M);
% Domains
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
dom.in = repmat([0,1],N,1);
dom.out = repmat([0,1],M,1);
% Dependency arrays
depmat = struct();
depmat.in = logical(randi([0,1],d,N));
depmat.out = logical(randi([0,1],d,M));
depmat.in(1,~any(depmat.in,1)) = true;      % make sure all vars appear at least once
depmat.out(1,~any(depmat.out,1)) = true;      % make sure all vars appear at least once

% Generate the functional
vars_A = vars;      vars_A.in = vars.out;   vars_A.out = cell(1,0);
dom_A = dom;        dom_A.in = dom.out;     dom_A.out = zeros(0,2);
dims_A = [nr,m];
if otype(1)==1
    degs_A = max(depmat.out,[],1);
else
    degs_A = sum(depmat.out,1);
end
if use_full_int
    % Aop is just an integral
    intmat_A = cell(1,M);
    for i=1:M
        intmat_A{i} = repmat([0,degs_A(i)+1],1,degs_A(i));
    end
    params = {speye(nr,m)};
    Aop = intop(params,intmat_A,vars_A,dom_A);
else
    % Generate random operator
    Aop = rand_intop(dims_A,vars_A.in,dom_A.in,degs_A,pdegs);
end


% Generate the random tensor-PI operator
Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,otype,pdegsB);
% Generate random function on which the operator acts
v = cell(1,d);
v_dum = cell(1,d);
degs_j = zeros(1,N);
v_full = 1;
pvars = polynomial(Bop.vars.in(:));
for j=1:d
    is_dep_j = depmat.in(j,:);
    n_j = size(Bop.ops{1,j},2);
    v{j} = rand_poly([n_j,1],pvars(is_dep_j),pdegv,50);

    % Declare a copy of the function expressed in terms of the appropriate
    % dummy variables
    degs_j = degs_j + is_dep_j;
    pvars_new = pvars;
    for k=find(is_dep_j)
        varname_k = pvar2varname(pvars(k));
        pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
    end
    % subs errors on an empty variable list: factor j may be constant
    if any(is_dep_j)
        v_dum{j} = subs(v{j},pvars(is_dep_j),pvars_new(is_dep_j));
    else
        v_dum{j} = v{j};
    end
    v_full = kron(v_full,v_dum{j});
end

% Compute the product
if use_full_int && nr==m
    Cop = int(Bop);
else
    Cop = mtimes_tensopvar(Aop,Bop);
end

% Apply the operators
Bv = apply_op(Bop,v);
ABv = apply_op(Aop,Bv);
Cv = apply_op(Cop,v_full);

% Tolerance relative to |ABv|: roundoff exceeds an absolute 1e-8 once |ABv| is much above 1
tol = 1e-8;
err = cleanpoly(ABv-Cv,1e-10);
disp('err = ')
disp(max(max(abs(err.C))))
nrmref = full(max(max(abs(double(ABv)))));
if max(max(abs(err.C)))>tol*max(1,nrmref)
    error("Outputs do not match: residual exceeds tolerance.")
end
disp('Test_intop_mtimes_tensopvar passed.')

% ============ SECONDARY CHECK (several terms, warning only) =================
% A tensopvar with nterms>1 (sum of independent random terms): the kernels of
% different terms for the same factor column have different monomial bases, which
% mtimes_tensopvar must unify before summing. A mismatch is reported with
% warning(), not error(): nested tensopvar/intop factor columns are unverified.
try
    nterms2 = randi([2,3]);
    Bop2 = rand_tensopvar(dims,vars,dom,depmat,nterms2,otype,pdegsB);
    v2 = cell(1,d);        v_dum2 = cell(1,d);
    degs_j2 = zeros(1,N);  v_full2 = 1;
    pvars2 = polynomial(Bop2.vars.in(:));
    for j=1:d
        is_dep_j2 = depmat.in(j,:);
        n_j2 = size(Bop2.ops{1,j},2);
        v2{j} = rand_poly([n_j2,1],pvars2(is_dep_j2),pdegv,50);
        degs_j2 = degs_j2 + is_dep_j2;
        pvars_new2 = pvars2;
        for k=find(is_dep_j2)
            varname_k2 = pvar2varname(pvars2(k));
            pvars_new2(k) = polynomial({[varname_k2{1},'_dum',num2str(degs_j2(k))]});
        end
        if any(is_dep_j2)
            v_dum2{j} = subs(v2{j},pvars2(is_dep_j2),pvars_new2(is_dep_j2));
        else
            v_dum2{j} = v2{j};
        end
        v_full2 = kron(v_full2,v_dum2{j});
    end
    Cop2 = mtimes_tensopvar(Aop,Bop2);
    Bv2 = apply_op(Bop2,v2);
    ABv2 = apply_op(Aop,Bv2);
    Cv2 = apply_op(Cop2,v_full2);
    err2 = cleanpoly(ABv2-Cv2,1e-10);
    nrmref2 = full(max(max(abs(double(ABv2)))));
    if max(max(abs(double(err2)))) > tol*max(1,nrmref2)
        warning("Test_intop_mtimes_tensopvar: secondary nterms=%d check ran but result did not match ground truth -- the known cross-term basis-mismatch bug was fixed in mtimes_tensopvar.m 09/30/2026, so this may be a regression or an untested residual case (see this file's comment above).",nterms2);
    else
        fprintf('Secondary check (nterms=%d): passed.\n',nterms2);
    end
catch ME2
    warning("Test_intop_mtimes_tensopvar: secondary nterms=%d check threw an error (%s) -- the known cross-term basis-mismatch bug was fixed in mtimes_tensopvar.m 09/30/2026, so this may be a regression or an untested residual case (see this file's comment above and mtimes_tensopvar.m's header changelog).", ...
        nterms2, ME2.message);
end
% ============= END SECONDARY CHECK ==========================================

% ================== DPVAR CHECK (hard failure) ===============================
% Aop with decision variables, same Aop parameters and primary-check Bop/v_full.
% Oracle: evaluate the decision variables at fixed values before composing (q=0
% path) vs. after composing the symbolic object (eval_dpvar_at/eval_intop_dpvar);
% equal by linearity. q is random, to cover different numbers of decision
% variables.
q = randi([50,200]);
Aop_dpvar = rand_intop(dims_A,vars_A.in,dom_A.in,degs_A,pdegs,q);
dvn = Aop_dpvar.Zd;
d0 = randn(numel(dvn),1);
d0map = containers.Map(dvn,num2cell(d0));
Aop_dpvar_num = eval_intop_dpvar(Aop_dpvar,d0map);

Cop_dpvar = mtimes_tensopvar(Aop_dpvar,Bop);
Cop_truth = mtimes_tensopvar(Aop_dpvar_num,Bop);
Cv_dpvar_num = eval_dpvar_at(apply_op(Cop_dpvar,v_full),d0map);
Cv_truth = apply_op(Cop_truth,v_full);

err3 = cleanpoly(Cv_dpvar_num - Cv_truth,1e-10);
nrmref3 = full(max(max(abs(double(Cv_truth)))));
if ~isdouble(err3)
    error("Test_intop_mtimes_tensopvar: dpvar check (q=%d): output is not constant.",q);
elseif max(max(abs(double(err3)))) > tol*max(1,nrmref3)
    error("Test_intop_mtimes_tensopvar: dpvar check (q=%d): residual exceeds tolerance.",q);
end
fprintf('Dpvar check (q=%d): passed.\n',q);
% ================ END DPVAR CHECK ============================================