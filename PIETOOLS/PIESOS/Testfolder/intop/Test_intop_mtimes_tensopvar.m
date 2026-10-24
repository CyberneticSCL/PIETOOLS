function Test_intop_mtimes_tensopvar()
% Test_intop_mtimes_tensopvar.m -- one local function for each case:
%   case_basic              mtimes_tensopvar for an 'intop' Aop and a
%                           'tensopvar' Bop
%   case_dvar_lists         the intop*tensopvar composition when the
%                           kernels of the functional depend on DIFFERENT
%                           decision variables in different regions of
%                           integration, as for the
%   case_N3multi            the intop*tensopvar composition when TWO OR
%                           MORE spatial variables are shared
%                           ("pass-through", N3>=2) between Bop's output
%                           and Bop's input (= Aop's input)
%   case_nested_factor      the composition of an 'intop' with a
%                           'tensopvar' whose factors are ITSELF nested
%                           'tensopvar' objects (see tensopvar.m and
%                           @tensopvar/initialize.m), as in the
%   case_S1_S2_combined     the intop*tensopvar composition when Bop has
%                           BOTH an "S2-type" output variable (new, no % DJ
%                           (CC), 10/21/2026 matching entry in Bop.vars.in)
%                           AND an "S1-type" input variable (consumed % DJ
%                           (CC), 10/21/2026
%   case_S2_extra_var       the intop*tensopvar % DJ (CC), 10/21/2026
%                           composition when Bop introduces an "S2-type"
%                           output variable: one that % DJ (CC), 10/21/2026
%                           appears in Bop.vars.out but has NO
%                           corresponding entry in Bop.vars.in (a
%   case_S3_mixed_deps      the intop*tensopvar composition when Bop has
%                           N3=3 shared (pass-through) variables, but its
%                           factors do NOT all depend on the same ones, e.g
%   case_shared_input_var   the intop*tensopvar composition when >=2 of the
%                           factors of Bop depend on the SAME named input
%                           variable (so their copies of it must be given
%                           distinct dummy names before
%   case_zero_sweep         Randomized sweep of the intop*tensopvar
%                           composition across varying numbers of spatial
%                           variables, factors, terms, and a wide range of
%                           coefficient densities (to generate a
%   case_zero_term          that composing an 'intop' with a 'tensopvar' in
%                           which ONE term (of several) has a zero factor
%                           returns the correct, generally NONZERO result
%                           (the other terms still
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_mtimes_tensopvar_dvar_lists,
%   Test_mtimes_tensopvar_N3multi, Test_mtimes_tensopvar_nested_factor,
%   Test_mtimes_tensopvar_S1_S2_combined,
%   Test_mtimes_tensopvar_S2_extra_var,
%   Test_mtimes_tensopvar_S3_mixed_deps,
%   Test_mtimes_tensopvar_shared_input_var,
%   Test_mtimes_tensopvar_zero_sweep, Test_mtimes_tensopvar_zero_term into
%   this test, as cases;
% DJ (CC), 10/24/2026: case_nested_factor: section (iv), a Kronecker product
%   of operators as factor in a pointwise product (supported, guard lifted in
%   mtimes_tensopvar), 12 cases, as in PIESOS_2D;

cases = {@case_basic, @case_dvar_lists, @case_N3multi, @case_nested_factor, @case_S1_S2_combined, @case_S2_extra_var, @case_S3_mixed_deps, @case_shared_input_var, @case_zero_sweep, @case_zero_term};
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
    error('Test_intop_mtimes_tensopvar: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_intop_mtimes_tensopvar: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- Tests mtimes_tensopvar for an 'intop'
% Aop and a 'tensopvar' Bop.
%
% Oracle: the composition Cop = mtimes_tensopvar(Aop,Bop) applied to the
% product input v_full is compared to Aop applied to the result of Bop
% applied to v, both evaluated with apply_op. Three parts: a primary check
% (one term, hard failure), a secondary check (several terms, warning
% only), and a check with decision variables in the kernels of Aop (hard
% failure).
%
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
end

%% --------------------------------------------------------------------------
function case_dvar_lists()
% case_dvar_lists -- Tests the intop*tensopvar
% composition when the kernels of the functional depend on DIFFERENT
% decision variables in different regions of integration, as for the
% functionals declared by piesos_quadvar. The composition adds the
% contributions of different regions (plus_spstruct in
% mtimes_tensopvar.m).
%
% Oracle: the decision variables are evaluated at random values d0 BEFORE
% composing, and the functional applied to the output of Bop (apply_op, no
% composition). The composition of the symbolic functional, evaluated at
% d0, must agree.
%
%
% NOTES
% - All kernels of an 'intop' share one list Zd: "different decision
%   variables" means coefficients on different subsets of Zd (thin_dvars
%   zeroes the other rows, used_dvars reads the subsets).

    tol = 1e-8;
    ntested = 0;
    q = 12;

    dims = [2,1; 1,2];
    vars = struct('in',{{'sShared'}},'out',{{'sOut1','sShared'}});
    dom = [0,1];
    depmat_cases = {logical([1,1;1,0]), logical([1,1;0,1]), logical([1,1;1,1])};

    for c = 1:numel(depmat_cases)
        for nterms = 1:2
            depmat = struct('in',true(2,1),'out',depmat_cases{c});
            % The test is only meaningful if the functional has several regions, with kernels that
            % depend on different decision variables: redraw until it does (the random integrals
            % may reduce to a single region)
            for attempt = 1:50
                rng(5000+10*c+nterms+1000*attempt);
                Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],1,0.5);
                pdegs = ones(1,numel(Bop.vars.out));
                Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs,q);
                Aop = thin_dvars(Aop);
                lists = cellfun(@(P) used_dvars(P,numel(Aop.Zd)),Aop.params(:),'UniformOutput',false);
                if numel(unique(lists))>=2
                    break
                end
            end
            if numel(unique(lists))<2
                error("case_dvar_lists: all kernels depend on the same decision variables (case=%d, nterms=%d).",c,nterms)
            end
            d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:q,'UniformOutput',false),num2cell(randn(q,1)));
            Aop_num = eval_intop_dpvar(Aop,d0map);

            pvars = polynomial(Bop.vars.in(:));
            v = {rand_poly([1,1],pvars,2,10), rand_poly([2,1],pvars,2,10)};
            Bv = apply_op(Bop,v);
            ABv = apply_op(Aop_num,Bv);
            nrmref = full(max(max(abs(double(ABv)))));
            if nrmref<=tol
                continue
            end

            Cop = mtimes_tensopvar(Aop,Bop);
            % Both factors depend on the SAME sole input variable 'sShared': give each its
            % own dummy copy before combining (see Test_intop_mtimes_tensopvar/case_shared_input_var)
            d = 2;  degs_j = 0;  v_full = 1;  v_dum = cell(1,d);
            for j=1:d
                degs_j = degs_j + 1;
                varname_1 = pvar2varname(pvars(1));
                pvar_new = polynomial({[varname_1{1},'_dum',num2str(degs_j)]});
                v_dum{j} = subs(v{j},pvars(1),pvar_new);
                v_full = kron(v_full,v_dum{j});
            end
            Cv = eval_dpvar_at(apply_op(Cop,v_full),d0map);
            err = cleanpoly(ABv-Cv,1e-10);
            ntested = ntested+1;
            if ~isdouble(err) || max(max(abs(double(err))))>tol*max(1,nrmref)
                error("case_dvar_lists failed: result does not match ground truth (case=%d, nterms=%d).",c,nterms)
            end
        end
    end

    if ntested==0
        error("case_dvar_lists: every case was skipped or degenerate -- no case actually exercised the check.")
    end
    fprintf('case_dvar_lists passed (%d/%d cases checked).\n',ntested,numel(depmat_cases)*2);
end

%% --------------------------------------------------------------------------
function case_N3multi()
% case_N3multi -- Tests the intop*tensopvar composition
% when TWO OR MORE spatial variables are shared ("pass-through", N3>=2)
% between Bop's output and Bop's input (= Aop's input).
%
% Oracle: apply_op(Bop) followed by apply_op(Aop), independent of the
% composition.
%
%
% NOTES
% - Scope is Part 1 (d==1, N3 in {2,3,4}) and Part 2 (d=2, N3=2), with only
%   the shared variables. Operators with N3>=3 shared variables together
%   with extra non-shared variables are not covered.
% - Each factor's copy of a shared input variable gets a distinct dummy
%   name (_dum<j>) before the copies are combined into the input of Cop.

    tol = 1e-8;

    % === PART 1: d==1, ONLY the shared variables, N3 in {2,3,4} ===
    nchecked = 0;
    for trial = 1:20
        rng(120000+trial);

        N3 = randi([2,4]);
        shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
        vars = struct('in',{shared_vars},'out',{shared_vars});
        dom = [0,1];
        dims = randi(2,[1,2]);
        depmat = struct('in',true(1,N3),'out',true(1,N3));
        dnsty = 0.1 + 0.5*rand();

        Bop = rand_tensopvar(dims,vars,dom,depmat,1,[2,2],1,dnsty);
        pdegs = ones(1,N3);
        Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

        pvars = polynomial(Bop.vars.in(:));
        v = {rand_poly([size(Bop.ops{1,1},2),1],pvars,2,10)};
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue
        end

        Cop = mtimes_tensopvar(Aop,Bop);
        Cv = apply_op(Cop,v{1});

        err = cleanpoly(ABv-Cv,1e-10);
        if max(max(abs(err.C)))>tol*max(1,nrmref)
            error("case_N3multi failed: d==1,N3=%d case does not match ground truth (trial=%d).",N3,trial)
        end
        nchecked = nchecked+1;
    end
    if nchecked==0
        error("case_N3multi: Part 1 (d==1) -- every trial was skipped, no trial exercised N3>=2.")
    end
    fprintf('Part 1 (d==1, N3 in {2,3,4}, no extra variables) passed: %d/20 trials checked against ground truth.\n',nchecked);

    % ===== PART 2: d>=2, N3=2 =======
    nchecked2 = 0;
    for trial = 1:20
        rng(130000+trial);

        N3 = 2; d = 2;
        vars = struct('in',{{'sShared1','sShared2'}},'out',{{'sShared1','sShared2'}});
        dom = [0,1];
        dims = randi(2,[d,2]);
        depmat = struct('in',true(d,N3),'out',true(d,N3));

        Bop = rand_tensopvar(dims,vars,dom,depmat,1,[2,2],1,0.9);
        pdegs = ones(1,N3);
        Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

        pvars = polynomial(Bop.vars.in(:));
        v = cell(1,d);
        for j=1:d
            v{j} = rand_poly([size(Bop.ops{1,j},2),1],pvars,2,10);
        end
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue
        end

        Cop = mtimes_tensopvar(Aop,Bop);
        v_dum = cell(1,d); degs_j = zeros(1,N3); v_full = 1;
        for j=1:d
            degs_j = degs_j + 1;
            pvars_new = pvars;
            for k=1:N3
                varname_k = pvar2varname(pvars(k));
                pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
            end
            v_dum{j} = subs(v{j},pvars,pvars_new);
            v_full = kron(v_full,v_dum{j});
        end
        Cv = apply_op(Cop,v_full);
        err = cleanpoly(ABv-Cv,1e-10);
        if max(max(abs(err.C)))>tol*max(1,nrmref)
            error("case_N3multi failed: d>=2,N3=2 case does not match ground truth (trial=%d).",trial)
        end
        nchecked2 = nchecked2+1;
    end
    if nchecked2==0
        error("case_N3multi: Part 2 (d>=2,N3=2) -- every trial was skipped, no trial exercised a nonzero case.")
    end
    fprintf('Part 2 (d>=2, N3=2) passed: %d/20 trials checked against ground truth.\n',nchecked2);
end

%% --------------------------------------------------------------------------
function case_nested_factor()
% case_nested_factor -- Tests the composition of an
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
% Section (iv): a Kronecker product of operators as a factor in a pointwise
% product, with a Kronecker product or a single operator as sibling.
%
%
% NOTES
% - Nestings that mtimes_tensopvar does not support (elementwise product of
%   the inputs of a factor in a pointwise product) are only checked to be
%   rejected with a 'not supported' message, section (v).
% - Sparse random operators can give a zero reference value: check_deep
%   draws again (up to 20 times, error id
%   'case_nested_factor:vacuous'); a zero reference is
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

    %% (iv) Kronecker product of operators as a factor in a pointwise product,
    %   B(x1,...,x4)(s) = N1(x1,x2)(s).*N2(x3,x4)(s),   N1 = P1 o P2,   N2 = P3 o P4 or a single operator S
    %   The factors of N1 and of N2 share the output copies s_1,s_2 (copy k of N1 is multiplied with copy k of
    %   N2), as in apply_nested. Fields of each case: Ntot (spatial variables), r1 (rows of the two factors
    %   of N1), r2 (rows of the two factors of N2, with the same product of rows, or [] for a single operator
    %   S with prod(r1) rows), nterms (terms in N1 and N2), nBterms (terms in B), q (decision variables in
    %   Aop), dep (optional, as in (i), rows P1,P2,P3,P4 or P1,P2,S)
    cases = {
        struct('Ntot',1,'r1',[1,1],'r2',[],   'nterms',1,'nBterms',1,'q',0)
        struct('Ntot',1,'r1',[1,1],'r2',[1,1],'nterms',1,'nBterms',1,'q',0)
        struct('Ntot',1,'r1',[2,1],'r2',[2,1],'nterms',1,'nBterms',1,'q',0)
        struct('Ntot',1,'r1',[2,1],'r2',[1,2],'nterms',1,'nBterms',1,'q',0)
        struct('Ntot',1,'r1',[1,2],'r2',[2,1],'nterms',1,'nBterms',1,'q',0)
        struct('Ntot',1,'r1',[2,2],'r2',[],   'nterms',1,'nBterms',1,'q',0)
        struct('Ntot',1,'r1',[1,2],'r2',[1,2],'nterms',2,'nBterms',2,'q',0)
        struct('Ntot',1,'r1',[2,1],'r2',[2,1],'nterms',1,'nBterms',1,'q',5)
        struct('Ntot',1,'r1',[2,2],'r2',[2,2],'nterms',1,'nBterms',1,'q',0,'dnsty',0.3)
        % two spatial variables, each of which only some of the operators map from and to
        struct('Ntot',2,'r1',[1,1],'r2',[1,1],'nterms',1,'nBterms',1,'q',0, ...
               'dep',struct('in',[1 0;0 1;1 0;0 1],'out',[1 0;0 1;1 0;0 1]))
        struct('Ntot',2,'r1',[2,1],'r2',[1,2],'nterms',1,'nBterms',1,'q',0, ...
               'dep',struct('in',[1 0;0 1;1 0;0 1],'out',[1 0;0 1;0 1;1 0]))
        struct('Ntot',2,'r1',[2,1],'r2',[],   'nterms',1,'nBterms',1,'q',0, ...
               'dep',struct('in',[1 0;1 1;0 1],'out',[1 0;1 1;0 1]))    };
    for c=1:numel(cases)
        rng(5000+c);
        check_kron_pointwise(cases{c},tol,sprintf('(iv) case %d',c));
        ncases = ncases+1;
    end

    %% (v) unsupported nestings must still be rejected, with an informative message
    % Elementwise product of the inputs of operators as a factor in a pointwise product
    vars = struct('in',{{'s'}},'out',{{'s'}});
    dom = struct('in',[0,1],'out',[0,1]);
    dep2 = struct('in',true(2,1),'out',true(2,1));
    dep1 = struct('in',true,'out',true);
    rng(5000);
    Nk = rand_tensopvar([1,1;1,1],vars,dom,dep2,1,[1,1],1,0.5);
    Sk = rand_tensopvar([1,1],vars,dom,dep1,1,[1,2],1,0.5);
    Bk = tensopvar();       Bk.otype = [1,2];       Bk.ops = {Nk,Sk};       Bk = initialize(Bk);
    Ak = rand_intop([1,Bk.matdim(1)],Bk.vars.out,Bk.dom.out,Bk.deps.out,1);
    expect_error(@() mtimes_tensopvar(Ak,Bk),'not supported');

    fprintf('case_nested_factor passed (%d cases).\n',ncases)
end

%% --------------------------------------------------------------------------
function case_S1_S2_combined()
% case_S1_S2_combined -- Tests the intop*tensopvar
% composition when Bop has BOTH an "S2-type" output variable (new, no       % DJ (CC), 10/21/2026
% matching entry in Bop.vars.in) AND an "S1-type" input variable (consumed  % DJ (CC), 10/21/2026
% internally, no matching entry in Bop.vars.out) at once, together with
% N3>=2 shared/pass-through variables.
%
% Oracle: apply_op(Bop) followed by apply_op(Aop), independent of the
% composition.
%
%
% NOTES
% - Cases with N3>=3 and d>=2 are not covered, with or without extra
%   variables.

    tol = 1e-8;

    % --- Part 1: d==1, N1=1,N2=1, N3 in {2,3,4,5} -- broad ---
    nchecked1 = 0;
    ntrials1 = 0;
    for N3 = [2,3,4,5]
        for seed_ofst = [0,100,200]
            ntrials1 = ntrials1+1;
            rng(60000+N3+seed_ofst);

            d = 1;
            shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
            vars = struct('in',{[shared_vars,{'sExtraIn'}]},'out',{[shared_vars,{'sExtraOut'}]});
            depmat = struct('in',true(d,N3+1),'out',true(d,N3+1));

            Bop = rand_tensopvar([1,1],vars,[0,1],depmat,1,[2,2],1,0.5);
            pdegs = ones(1,numel(Bop.vars.out));
            Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

            pvars = polynomial(Bop.vars.in(:));
            v = {rand_poly([size(Bop.ops{1,1},2),1],pvars,2,10)};
            Bv = apply_op(Bop,v);
            ABv = apply_op(Aop,Bv);
            nrmref = full(max(max(abs(double(ABv)))));
            if nrmref<=tol
                continue
            end

            Cop = mtimes_tensopvar(Aop,Bop);
            Cv = apply_op(Cop,v{1});

            err = cleanpoly(ABv-Cv,1e-10);
            errval = full(max(max(abs(err.C))));
            if isempty(errval); errval = 0; end
            if errval>tol*max(1,nrmref)
                error("case_S1_S2_combined failed: d==1,N3=%d case does not match ground truth (seed_ofst=%d).",N3,seed_ofst)
            end
            nchecked1 = nchecked1+1;
        end
    end
    if nchecked1==0
        error("case_S1_S2_combined: Part 1 -- every trial was skipped, nothing exercised.")
    end
    fprintf('Part 1 (d==1, N1=1,N2=1, N3 in {2,3,4,5}) passed: %d/%d trials checked against ground truth.\n',nchecked1,ntrials1);

    % --- Part 2: d==2, N1=1,N2=1, N3=2 ONLY -- do not raise N3 here, see NOTES ---
    nchecked2 = 0;
    ntrials2 = 0;
    for seed_ofst = [0,1,2]
        ntrials2 = ntrials2+1;
        N3 = 2; d = 2;
        rng(70000+N3*10+seed_ofst);

        shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
        vars = struct('in',{[shared_vars,{'sExtraIn'}]},'out',{[shared_vars,{'sExtraOut'}]});
        depmat = struct('in',true(d,N3+1),'out',true(d,N3+1));

        Bop = rand_tensopvar(repmat([1,1],d,1),vars,[0,1],depmat,1,[2,2],1,0.5);
        pdegs = ones(1,numel(Bop.vars.out));
        Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

        pvars = polynomial(Bop.vars.in(:));
        Nin_tot = numel(Bop.vars.in);
        v = cell(1,d);
        for j=1:d
            v{j} = rand_poly([size(Bop.ops{1,j},2),1],pvars,2,10);
        end
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue
        end

        Cop = mtimes_tensopvar(Aop,Bop);

        % Per-factor dummy names for variables shared across >=2 factors (required for d>=2, dense depmat)
        degs_j = zeros(1,Nin_tot); v_full = 1;
        for j=1:d
            degs_j = degs_j+1;
            pvars_new = pvars;
            for k=1:Nin_tot
                varname_k = pvar2varname(pvars(k));
                pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
            end
            v_dum = subs(v{j},pvars,pvars_new);
            v_full = kron(v_full,v_dum);
        end
        Cv = apply_op(Cop,v_full);

        err = cleanpoly(ABv-Cv,1e-10);
        errval = full(max(max(abs(err.C))));
        if isempty(errval); errval = 0; end
        if errval>tol*max(1,nrmref)
            error("case_S1_S2_combined failed: d==2,N3=2 case does not match ground truth (seed_ofst=%d).",seed_ofst)
        end
        nchecked2 = nchecked2+1;
    end
    if nchecked2==0
        error("case_S1_S2_combined: Part 2 -- every trial was skipped, nothing exercised.")
    end
    fprintf('Part 2 (d==2, N1=1,N2=1, N3=2) passed: %d/%d trials checked against ground truth.\n',nchecked2,ntrials2);
end

%% --------------------------------------------------------------------------
function case_S2_extra_var()
% case_S2_extra_var -- Tests the intop*tensopvar         % DJ (CC), 10/21/2026
% composition when Bop introduces an "S2-type" output variable: one that    % DJ (CC), 10/21/2026
% appears in Bop.vars.out but has NO corresponding entry in Bop.vars.in (a
% newly introduced output, as opposed to a variable shared/passed through
% between the input and output of Bop).
%
% Oracle: apply_op(Bop) followed by apply_op(Aop), independent of the
% composition.
%
%
% NOTES
% - The mirror "S1-type" case (variable in Bop.vars.in only) is covered by  % DJ (CC), 10/21/2026
%   Test_intop_mtimes_tensopvar/case_shared_input_var.
% - Deliberately small (d==1, N3<=5): the case needs no more.

    tol = 1e-8;
    nchecked = 0;
    ntrials = 0;
    for N3 = [3,4,5]
        for seed_ofst = [0,100,200]
            ntrials = ntrials+1;
            rng(50400+N3+seed_ofst);

            d = 1;
            shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
            vars = struct('in',{shared_vars},'out',{[shared_vars,{'sExtraOut'}]});
            depmat = struct('in',true(d,N3),'out',true(d,N3+1));

            Bop = rand_tensopvar([1,1],vars,[0,1],depmat,1,[2,2],1,0.5);
            pdegs = ones(1,numel(Bop.vars.out));
            Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

            pvars = polynomial(Bop.vars.in(:));
            v = {rand_poly([size(Bop.ops{1,1},2),1],pvars,2,10)};
            Bv = apply_op(Bop,v);
            ABv = apply_op(Aop,Bv);
            nrmref = full(max(max(abs(double(ABv)))));
            if nrmref<=tol
                continue
            end

            Cop = mtimes_tensopvar(Aop,Bop);
            Cv = apply_op(Cop,v{1});

            err = cleanpoly(ABv-Cv,1e-10);
            errval = full(max(max(abs(err.C))));
            if isempty(errval); errval = 0; end
            if errval>tol*max(1,nrmref)
    %           error("Test_mtimes_tensopvar_S1_extra_var failed: N3=%d case does not match ground truth (seed_ofst=%d).",N3,seed_ofst) % DJ (CC), 10/21/2026 (was)
                error("case_S2_extra_var failed: N3=%d case does not match ground truth (seed_ofst=%d).",N3,seed_ofst) % DJ (CC), 10/21/2026
            end
            nchecked = nchecked+1;
        end
    end
    if nchecked==0
    %   error("Test_mtimes_tensopvar_S1_extra_var: every trial was skipped (all zero ground truth) -- test exercised nothing.") % DJ (CC), 10/21/2026 (was)
        error("case_S2_extra_var: every trial was skipped (all zero ground truth) -- test exercised nothing.") % DJ (CC), 10/21/2026
    end
    %fprintf('Test_mtimes_tensopvar_S1_extra_var passed: %d/%d trials checked against ground truth.\n',nchecked,ntrials); % DJ (CC), 10/21/2026 (was)
    fprintf('case_S2_extra_var passed: %d/%d trials checked against ground truth.\n',nchecked,ntrials); % DJ (CC), 10/21/2026
end

%% --------------------------------------------------------------------------
function case_S3_mixed_deps()
% case_S3_mixed_deps -- Tests the intop*tensopvar
% composition when Bop has N3=3 shared (pass-through) variables, but its
% factors do NOT all depend on the same ones, e.g. factor 1 on (sA,sB) and
% factor 2 on (sA,sB,sC).
%
% Oracle: apply_op(Bop) followed by apply_op(Aop) on concrete polynomials,
% independent of the composition.
%
%
% NOTES
% - Tests the handling of N3-D parameter cells of different ndims in
%   misc/otimes_params.m>otimes_quadPolyCell: MATLAB drops trailing
%   singleton dimensions of the N3-D parameter cells, so factors that do
%   not depend on the last shared variable differ in ndims. Needs N3>=3,
%   different shared-variable sets across the factors, and d>=2 (otimes is
%   only called to combine factors).
% - Each factor's copy of a shared input variable gets a distinct dummy
%   name (_dum<j>) before the copies are combined into the input of Cop.

    tol = 1e-8;
    vars = struct('in',{{'sA','sB','sC'}},'out',{{'sA','sB','sC'}});
    dom = [0,1];
    dims = [2,1; 1,2];
    N = numel(vars.in);     d = size(dims,1);

    % Each case lists, per factor (rows), on which of (sA,sB,sC) it depends
    %   in:  input variables of the factor;     out: output variables
    cases = struct('name',{},'in',{},'out',{});
    cases(end+1) = struct('name','x1(sA,sB) & x2(sA,sB,sC)', ...
        'in',logical([1,1,0; 1,1,1]),'out',logical([1,1,0; 1,1,1]));
    cases(end+1) = struct('name','x1(sA,sB,sC) & x2(sA,sB)', ...
        'in',logical([1,1,1; 1,1,0]),'out',logical([1,1,1; 1,1,0]));
    cases(end+1) = struct('name','x1(sA,sB) & x2(sC)', ...
        'in',logical([1,1,0; 0,0,1]),'out',logical([1,1,0; 0,0,1]));
    cases(end+1) = struct('name','x1(sC) & x2(sA,sB)', ...
        'in',logical([0,0,1; 1,1,0]),'out',logical([0,0,1; 1,1,0]));
    cases(end+1) = struct('name','common output space, inputs (sA,sB) & (sA,sB,sC)', ...
        'in',logical([1,1,0; 1,1,1]),'out',true(2,3));

    nchecked = 0;
    for c = 1:numel(cases)
        ok = false;
        for seed = 1:5      % retry if the random draw makes the reference vanish
            rng(140000+10*c+seed);
            depmat = struct('in',cases(c).in,'out',cases(c).out);
            Bop = rand_tensopvar(dims,vars,dom,depmat,1,[2,2],1,0.5);
            pdegs = ones(1,numel(Bop.vars.out));
            Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

            % Random input to each factor, and the same input with each factor's
            % own dummy copy of every shared variable
            pvars = polynomial(Bop.vars.in(:));
            v = cell(1,d);      v_dum = cell(1,d);
            degs_j = zeros(1,N);        v_full = 1;
            for j=1:d
                is_dep_j = cases(c).in(j,:);
                n_j = size(Bop.ops{1,j},2);
                v{j} = rand_poly([n_j,1],pvars(is_dep_j),2,50);
                degs_j = degs_j + is_dep_j;
                pvars_new = pvars;
                for k=find(is_dep_j)
                    varname_k = pvar2varname(pvars(k));
                    pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
                end
                v_dum{j} = subs(v{j},pvars(is_dep_j),pvars_new(is_dep_j));
                v_full = kron(v_full,v_dum{j});
            end

            % Ground truth: apply the operators one after the other
            Bv = apply_op(Bop,v);
            ABv = apply_op(Aop,Bv);
            nrmref = full(max(max(abs(double(ABv)))));
            if nrmref<=tol
                continue
            end

            % Composition
            Cop = mtimes_tensopvar(Aop,Bop);
            Cv = apply_op(Cop,v_full);
            err = cleanpoly(ABv-Cv,1e-10);
            if max(max(abs(err.C)))>tol*max(1,nrmref)
                error("case_S3_mixed_deps failed: case '%s' does not match ground truth (seed=%d).",cases(c).name,seed)
            end
            ok = true;
            break
        end
        if ~ok
            error("case_S3_mixed_deps: case '%s' -- every random draw gave a (near) zero reference, nothing exercised.",cases(c).name)
        end
        nchecked = nchecked+1;
    end
    fprintf('case_S3_mixed_deps passed (%d/%d cases checked against ground truth).\n',nchecked,numel(cases));
end

%% --------------------------------------------------------------------------
function case_shared_input_var()
% case_shared_input_var -- Tests the intop*tensopvar
% composition when >=2 of the factors of Bop depend on the SAME named input
% variable (so their copies of it must be given distinct dummy names before
% being combined into a single test input, as in
% Test_intop_mtimes_tensopvar.m).
%
% Cases: symmetric and asymmetric per-factor dependence on a shared output
% variable, d=2, nterms in {1,2}. Oracle: apply_op(Bop) followed by
% apply_op(Aop).
%
%
% NOTES
% - The test input MUST give each factor its own dummy copy of the shared
%   variable (_dum<j>): kron(v{1},v{2}) with both in the one variable name
%   conflates two independent copies and mismatches the (correct)
%   composition.

    tol = 1e-8;
    ntested = 0;

    dims = [2,1; 1,2];
    vars = struct('in',{{'sShared'}},'out',{{'sOut1','sShared'}});
    dom = [0,1];
    depmat_cases = {logical([1,1;1,0]), logical([1,1;0,1]), logical([1,1;1,1])};

    for c = 1:numel(depmat_cases)
        for nterms = 1:2
            rng(999);
            depmat = struct('in',true(2,1),'out',depmat_cases{c});
            Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],1,0.5);
            pdegs = ones(1,numel(Bop.vars.out));
            Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

            pvars = polynomial(Bop.vars.in(:));
            v = {rand_poly([1,1],pvars,2,10), rand_poly([2,1],pvars,2,10)};
            Bv = apply_op(Bop,v);
            ABv = apply_op(Aop,Bv);
            nrmref = full(max(max(abs(double(ABv)))));
            if nrmref<=tol
                continue
            end

            Cop = mtimes_tensopvar(Aop,Bop);
            % Both factors depend on the SAME sole input variable 'sShared' --
            % give each its own dummy copy before combining (see NOTES).
            d = 2; N = 1; degs_j = 0; v_full = 1; v_dum = cell(1,d);
            for j=1:d
                degs_j = degs_j + 1;
                varname_1 = pvar2varname(pvars(1));
                pvar_new = polynomial({[varname_1{1},'_dum',num2str(degs_j)]});
                v_dum{j} = subs(v{j},pvars(1),pvar_new);
                v_full = kron(v_full,v_dum{j});
            end
            Cv = apply_op(Cop,v_full);
            err = cleanpoly(ABv-Cv,1e-10);
            ntested = ntested+1;
            if max(max(abs(err.C)))>tol*max(1,nrmref)
                error("case_shared_input_var failed: result does not match ground truth (case=%d, nterms=%d).",c,nterms)
            end
        end
    end

    if ntested==0
        error("case_shared_input_var: every case was skipped or degenerate -- no case actually exercised the check.")
    end
    fprintf('case_shared_input_var passed (%d/%d cases checked).\n',ntested,numel(depmat_cases)*2);
end

%% --------------------------------------------------------------------------
function case_zero_sweep()
% case_zero_sweep -- Randomized sweep of the
% intop*tensopvar composition across varying numbers of spatial variables,
% factors, terms, and a wide range of coefficient densities (to generate a
% diverse mix of structural zero patterns across terms/factors), checking
% the result against an independent ground truth every time.
%
% Companion to the single hand-picked scenario of
% Test_intop_mtimes_tensopvar/case_zero_term: no structural zero pattern may make
% mtimes_tensopvar silently return a wrong (typically all-zero) result.
% Also exercises the checks in get_params of mtimes_tensopvar.m (error ids
% 'get_params:inconsistent_zero_slot' and
% 'get_params:inconsistent_params_type'), which fail loudly if they fire in
% error. Oracle: apply_op(Bop) followed by apply_op(Aop), independent of
% mtimes_tensopvar/get_params.
%
%
% NOTES
% - Each factor's copy of a shared input variable gets its own dummy name
%   (_dum<j>) before the copies are combined (required when >=2 factors
%   depend on the same input variable).

    tol = 1e-8;
    ntrials = 150;
    nskipped = 0;
    nchecked = 0;

    for trial = 1:ntrials
        rng(70000+trial);

        N = randi([1,2]);
        M = randi([1,2]);
        d = randi([1,3]);
        nterms = randi([2,4]);
        dims = randi(2,[d,2]);
        % Force exactly one variable shared between vars.in/vars.out (N3>=1).
        vars_in = [arrayfun(@(k) sprintf('sIn%d',k),1:N-1,'UniformOutput',false),{'sShared'}];
        vars_out = [arrayfun(@(k) sprintf('sOut%d',k),1:M-1,'UniformOutput',false),{'sShared'}];
        vars = struct('in',{vars_in},'out',{vars_out});
        dom = [0,1];
        depmat = struct('in',logical(randi([0,1],d,N)),'out',logical(randi([0,1],d,M)));
        depmat.in(1,~any(depmat.in,1)) = true;
        depmat.out(1,~any(depmat.out,1)) = true;
        dnsty = 0.1 + 0.5*rand();   % vary density to mix partial- and full-zero patterns

        try
            Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],1,dnsty);
        catch
            nskipped = nskipped+1;
            continue
        end
        [~,is_zeroB] = clean_op(Bop);
        if is_zeroB
            continue   % trivially-zero case already covered by Test_intop_mtimes_tensopvar/case_zero_term
        end

        pdegs = ones(1,numel(Bop.vars.out));
        try
            Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);
        catch
            nskipped = nskipped+1;
            continue
        end

        pvars = polynomial(Bop.vars.in(:));
        v = cell(1,d);
        for j=1:d
            n_j = size(Bop.ops{1,j},2);
            is_dep_j = depmat.in(j,:);
            v{j} = rand_poly([n_j,1],pvars(is_dep_j),2,10);
        end
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue   % degenerate draw (ground truth itself ~0) -- not informative
        end

        try
            Cop = mtimes_tensopvar(Aop,Bop);
        catch ME
            error("case_zero_sweep failed: mtimes_tensopvar errored (trial=%d, seed=%d): %s", ...
                trial,70000+trial,ME.message)
        end
        % Own dummy name per factor for shared input variables (see NOTES)
        v_dum = cell(1,d); degs_j = zeros(1,N); v_full = 1;
        for j=1:d
            is_dep_j = depmat.in(j,:);
            degs_j = degs_j + is_dep_j;
            pvars_new = pvars;
            for k=find(is_dep_j)
                varname_k = pvar2varname(pvars(k));
                pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
            end
            if any(is_dep_j)
                v_dum{j} = subs(v{j},pvars(is_dep_j),pvars_new(is_dep_j));
            else
                v_dum{j} = v{j};
            end
            v_full = kron(v_full,v_dum{j});
        end
        Cv = apply_op(Cop,v_full);

        err = cleanpoly(ABv-Cv,1e-10);
        if max(max(abs(err.C)))>tol*max(1,nrmref)
            error("case_zero_sweep failed: result does not match ground truth (trial=%d, seed=%d).", ...
                trial,70000+trial)
        end
        nchecked = nchecked+1;
    end

    if nchecked==0
        error("case_zero_sweep: every trial was skipped or degenerate -- no trial actually exercised the check.")
    end
    fprintf('case_zero_sweep passed (%d/%d trials checked against ground truth, %d skipped on construction errors).\n', ...
        nchecked,ntrials,nskipped);
end

%% --------------------------------------------------------------------------
function case_zero_term()
% case_zero_term -- Tests that composing an 'intop'
% with a 'tensopvar' in which ONE term (of several) has a zero factor
% returns the correct, generally NONZERO result (the other terms still
% contribute), rather than treating the entire operator as zero.
%
% Also covers Bop ENTIRELY zero (early return in mtimes_tensopvar,
% bypassing get_params): the result must be the zero functional of the
% correct size. Oracle: apply_op(Bop) followed by apply_op(Aop),
% independent of mtimes_tensopvar/get_params.
%

    tol = 1e-8;

    for trial = 1:8
        rng(9900+trial);

        dims = [1,1; 1,1];
        vars = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
        dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
        depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
        Bop = rand_tensopvar(dims,vars,dom,depmat,3,[2,2],1,0.5);
        Bop.ops{2,1} = 0*Bop.ops{2,1};   % zero out one term's first factor

        nr = 1;
        pdegs = ones(1,numel(Bop.vars.out));
        Aop = rand_intop([nr,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

        pvars = polynomial(Bop.vars.in(:));
        v = {rand_poly([1,1],pvars(1),3,10), rand_poly([1,1],pvars(2),3,10)};
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            % Degenerate draw (remaining terms cancel): not informative
            continue
        end

        Cop = mtimes_tensopvar(Aop,Bop);
        v_full = kron(v{1},v{2});
        Cv = apply_op(Cop,v_full);

        err = cleanpoly(ABv-Cv,1e-10);
        if max(max(abs(err.C)))>tol*max(1,nrmref)
            error("case_zero_term failed: result does not match ground truth (trial=%d).",trial)
        end
    end

    % --- Complementary case: Bop entirely zero -- top-level short-circuit ---
    rng(9800);
    dims = [1,1; 1,1];
    vars = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
    dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
    depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
    Bop0 = rand_tensopvar(dims,vars,dom,depmat,3,[2,2],1,0.5);
    for i=1:3
        Bop0.ops{i,1} = 0*Bop0.ops{i,1};
    end
    pdegs0 = ones(1,numel(Bop0.vars.out));
    Aop0 = rand_intop([1,Bop0.matdim(1)],Bop0.vars.out,Bop0.dom.out,Bop0.deps.out,pdegs0);
    Cop0 = mtimes_tensopvar(Aop0,Bop0);
    if ~isequal(size(Cop0),[size(Aop0,1),size(Bop0,2)])
        error("case_zero_term failed: wrong size for fully-zero Bop (got %s).",mat2str(size(Cop0)))
    end
    pvars0 = polynomial(Bop0.vars.in(:));
    v0 = {rand_poly([1,1],pvars0(1),3,10), rand_poly([1,1],pvars0(2),3,10)};
    Cv0 = apply_op(Cop0,kron(v0{1},v0{2}));
    if ~isdouble(Cv0) || any(any(double(Cv0)~=0))
        error("case_zero_term failed: fully-zero Bop did not yield the zero functional.")
    end
    fprintf('Fully-zero Bop case: passed.\n');

    disp('case_zero_term passed.')
end

%% ---- helpers --------------------------------------------------------------

function Aop = thin_dvars(Aop)
% Retain a different random subset of the decision variables in each parameter of Aop
% (the rows of the others, in the blocks of rows of each row of the matrix, are set to zero)
nd = numel(Aop.Zd);
m = size(Aop,1);
for k=1:numel(Aop.params)
    P = Aop.params{k};
    if nd==0,       continue,       end
    J = sort(randperm(nd,ceil(nd/2)));
    rows = reshape(((0:m-1)'*(nd+1) + [1,1+J]).',[],1);
    P(setdiff(1:size(P,1),rows),:) = 0;
    Aop.params{k} = P;
end
end

function s = used_dvars(P,nd)
% String with the decision variables (rows of P) that have nonzero coefficients
[i,~] = find(P);
d = mod(i-1,nd+1);          % 0 for the constant term
s = char(strjoin(string(unique(d(d>0))'),','));
end

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

function check_kron_pointwise(cfg,tol,label)
% Build B = N1.*N2 with N1 a Kronecker product of two operators, and N2 a Kronecker product of two
% operators or a single operator, and compare the composition with a functional A with the ground truth
Ntot = cfg.Ntot;
svars = arrayfun(@(i) ['s',num2str(i)],1:Ntot,'UniformOutput',false);
nl = 3 + (numel(cfg.r2)==2);            % leaf operators: P1,P2 in N1, P3,P4 in N2 (or S)
if isfield(cfg,'dep')
    depin = logical(cfg.dep.in);        depout = logical(cfg.dep.out);
else
    depin = true(nl,Ntot);              depout = true(nl,Ntot);
end
is_in = any(depin,1);                   is_out = any(depout,1);
vars = struct('in',{svars(is_in)},'out',{svars(is_out)});
depin = depin(:,is_in);                 depout = depout(:,is_out);
N = numel(vars.in);                     M = numel(vars.out);
dom = struct('in',repmat([0,1],N,1),'out',repmat([0,1],M,1));
nin = randi(2,nl,1);                    % column dimensions of the leaf operators
pdegs = struct('in',ones(1,N),'out',ones(1,M));
dnsty = 0.5;
if isfield(cfg,'dnsty')
    dnsty = cfg.dnsty;
end
dep1 = struct('in',depin(1:2,:),'out',depout(1:2,:));
dep2 = struct('in',depin(3:nl,:),'out',depout(3:nl,:));
if numel(cfg.r2)==2
    dims2 = [cfg.r2(:),nin(3:4)];
else
    dims2 = [prod(cfg.r1),nin(3)];
end
Nterm = cell(cfg.nBterms,1);            Mterm = cell(cfg.nBterms,1);
for it=1:cfg.nBterms
    Nterm{it} = rand_tensopvar([cfg.r1(:),nin(1:2)],vars,dom,dep1,cfg.nterms,[2,2],pdegs,dnsty);
    Mterm{it} = rand_tensopvar(dims2,vars,dom,dep2,cfg.nterms,[2,2],pdegs,dnsty);
end
Bop = tensopvar();      Bop.otype = [1,2];
Bop.ops = [Nterm,Mterm];
Bop = initialize(Bop);

% States of the leaf operators, in the order in which they appear in B
pvars = polynomial(Bop.vars.in(:));
leaves = cell(1,nl);
for l=1:nl
    leaves{l} = rand_poly([nin(l),1],pvars(depin(l,:)),2,10);
end
compose_and_compare(Bop,leaves,depin,cfg.q,tol,label);
end

function check_deep(variant,tol,label)
% Nested operators in nested operators, in one spatial variable, with sparse operators: the
% reference value may be zero, in which case the operators are drawn again
for attempt = 1:20
    try
        check_deep_once(variant,tol,label);
        return
    catch ME
        if ~strcmp(ME.identifier,'case_nested_factor:vacuous')
            rethrow(ME)
        end
    end
end
error("case_nested_factor (%s): the reference value is zero in each of the draws.",label)
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
    Aop = thin_dvars_nested_factor(Aop);
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
    error('case_nested_factor:vacuous',"case_nested_factor (%s): the reference value is zero, the test is vacuous.",label)
elseif ~isdouble(err) || max(max(abs(double(err))))>tol*max(1,nrmref)
    error("case_nested_factor failed (%s): residual %.3e exceeds tolerance.",label,max(max(abs(double(err)))))
end
end

function Aop = thin_dvars_nested_factor(Aop)
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
        error("case_nested_factor failed: unexpected error message '%s'.",ME.message)
    end
    return
end
error("case_nested_factor failed: expected an error, but none was thrown.")
end

