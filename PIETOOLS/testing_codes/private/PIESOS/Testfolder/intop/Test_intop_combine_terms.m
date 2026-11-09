function Test_intop_combine_terms()
% Test_intop_combine_terms.m -- one local function for each case:
%   case_basic              @intop/combine_terms.m: reordering the dummy
%                           variables of a functional acting on a
%                           distributed monomial of nx states over N=2
%                           spatial variables must not change its action.
%                           Oracle:
%   case_2D_twisted         @intop/combine_terms.m (and
%                           @intop/reorder_states.m, @intop/reorder_vars.m,
%                           which it calls) with a "twisted" nesting order:
%                           the two spatial directions s1, s2 integrate
%   case_canonical          that @intop/combine_terms.m (the form with
%                           arguments (Kop,degs,dims,depmat)) gives the
%                           same intmat and the same kernels for equal
%                           functionals, so that the coefficients of
%   case_idempotent         that @intop/combine_terms.m (the form with
%                           arguments (Kop,degs,dims,depmat)) puts the
%                           repeated factors of a state in a canonical
%                           order, so that
%   case_ND_random          Randomized test of @intop/combine_terms.m for a
%                           repeated state spanning N>=2 spatial
%                           directions, complementing the hand-picked
%                           scenarios in
%   case_sum_blocks         that @intop/combine_terms.m, which sums the
%                           blocks of terms with the same order of the
%                           copies of a state in one pass, gives a
%                           functional that
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_intop_combine_terms_2D_twisted,
%   Test_intop_combine_terms_canonical,
%   Test_intop_combine_terms_idempotent,
%   Test_intop_combine_terms_ND_random, Test_intop_combine_terms_sum_blocks
%   into this test, as cases;

cases = {@case_basic, @case_2D_twisted, @case_canonical, @case_idempotent, @case_ND_random, @case_sum_blocks};
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
    error('Test_intop_combine_terms: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_intop_combine_terms: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- Tests @intop/combine_terms.m: reordering
% the dummy variables of a functional acting on a distributed monomial of
% nx states over N=2 spatial variables must not change its action. Oracle:
% apply_op of the functional before vs. after combine_terms, on the same
% random monomial.
%

% NOTES
% - nx=2 (degree-2 states, di=4 in each direction) is used because
%   rand_intop cannot generate the functional for nx=3 and depmat=true
%   (di=6 in each direction): split_ints then expands each term into up to
%   6!=720 sub-terms per direction, exceeding the maximal array size.
% - Value check only (relative residual); the canonical form of the result
%   is tested in Test_intop_combine_terms/case_canonical.

% Set parameters defining the distributed monomial

    N = 2;                      % number of spatial variables
    nx = 2;                      % number of state variables
    dmin = 2;
    dmax = 2;                   % maximal distributed monomial degree
    degs = randi([dmin,dmax],1,nx);    % degree of each state in the monomial
    d = sum(degs);
    state_idcs = repelem(1:nx,degs);
    n_arr = randi(2,[1,nx]);     % row dimension of state variables
    q_arr = ones(1,nx);          % column dimensions of state variables (should be 1)
    depmat = true(nx,N);
    n_arr_in = n_arr(state_idcs);
    q_arr_in = q_arr(state_idcs);
    depmat_in = depmat(state_idcs,:);

    % Set parameters defining functional
    m = 2;  
    n = prod(n_arr_in);            % column dimension of functional       
    dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
    deps = sum(depmat_in,1);
    pdegs = 3;

    % Set parameters defining the function
    q = 2;              % column dimension of function
    pdegv = 2;

    % Generate a random intop object
    Aop = rand_intop([m,n],N,dom,deps,pdegs);
    Aop2 = Aop;

    % Reorder the variables
    Anew = combine_terms(Aop2,degs,n_arr,depmat);

    % Generate a random distributed monomial on which the functional acts
    vnames = cell(1,N);
    for i=1:N
        vnames{i} = ['s',num2str(i)];
    end
    xcell = cell(1,nx);
    pvars_arr = cell(1,nx);
    for i=1:nx
        dep_i = depmat(i,:);
        vnames_i = vnames(dep_i);
        pvars_arr{i} = polynomial(vnames_i(:))';
        xcell{i} = rand_poly([n_arr(i),q_arr(i)],pvars_arr{i},pdegv,50);
        if all(isequal(xcell{i},0))
            xcell{i} = eye(n_arr(i),q_arr(i));
        end
    end
    vnames_in = Aop.varnames.in;
    ycell = cell(1,d);
    xval = 1;
    cntr = zeros(1,N);
    pvars_new = cell(1,d);
    for j=1:d
        dep_j = depmat_in(j,:);
        vnames_j = vnames_in(dep_j);
        cntr(dep_j) = cntr(dep_j)+1;    % move to the next variable
        vars_j = cellfun(@(a,b) a(b),vnames_j,num2cell(cntr(dep_j)),'UniformOutput',false);
        pvars_new{j} = polynomial([vars_j{:}])';
        ycell{j} = subs(xcell{state_idcs(j)},pvars_arr{state_idcs(j)},pvars_new{j});
        xval = kron(xval,ycell{j});
    end

    % Apply the functionals
    Ax = apply_op(Aop2,xval);
    Ax_alt = apply_op(Anew,xval);


    % Compare
    tol = 1e-8;
    err_A = cleanpoly(Ax - Ax_alt,1e-10);
    if ~isdouble(err_A)
        error("At least one of the outputs is not constant.")
    elseif norm(double(err_A),'fro')/norm(double(Ax),'fro')>tol
        error("Outputs do not match: relative residual exceeds tolerance.")
    else
        disp('err_C = ')
        disp(double(err_A)/norm(double(Ax)))
    end
end

%% --------------------------------------------------------------------------
function case_2D_twisted()
% case_2D_twisted -- Tests @intop/combine_terms.m
% (and @intop/reorder_states.m, @intop/reorder_vars.m, which it calls) with
% a "twisted" nesting order: the two spatial directions s1, s2 integrate
% the SAME pair of repeated factor-slots in OPPOSITE order (slot1 inner/
% slot2 outer along s1, but slot1 outer/slot2 inner along s2). Parts 4/5
% are multi-group cases (two distinct nesting orders along s1, so
% combine_terms splits into two groups, canonicalizes each with a
% different permutation, then sums). Oracle: apply_op of the functional
% before vs. after combine_terms on x(s1,s2) o x(s1,s2) (o x(s1,s2) for
% degree 3), exact polynomial arithmetic.
%

% NOTES
% - Every copy of the state must be the same base formula, renamed by
%   subs() to its own dummy variables (s1_1/s2_1, s1_2/s2_2, ...):
%   apply_op matches dummy variables to slots by NAME, not by kron()
%   position, and combine_terms assumes repeated slots hold the same
%   function (one state to a tensor power). Independent functions per slot
%   violate that precondition.
% - Value check only; the canonical form is tested in
%   Test_intop_combine_terms/case_canonical.

    tol = 1e-8;
    vars = struct('in',{{'s1','s2'}},'out',{cell(1,0)});
    dom = struct('in',[0,1;0,1],'out',zeros(0,2));
    depmat = true(1,2);   % the (single) state depends on both s1 and s2

    %% Part 1: single-group, SCALAR state, degree 2, twisted order
    I1 = [0,2, 0,3];   % s1: slot1 inner [0,dum2], slot2 outer [0,b]
    I2 = [0,3, 0,1];   % s2: slot1 outer [0,b], slot2 inner [0,dum1]
    intmat = {I1,I2};

    pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
    K = s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*s1_dum2^2*s2_dum1;
    Kop = intop({K},intmat,vars,dom,false);
    Kop2 = combine_terms(Kop,2,1,depmat);

    pvar s t s1_1 s2_1 s1_2 s2_2
    vx = s^2*t + 3*t^3 - s*t^2 + 1;
    xval = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);

    Y0 = apply_op(Kop,xval);
    Y1 = apply_op(Kop2,xval);
    err = cleanpoly(Y0-Y1,1e-10);
    if ~isdouble(err) || abs(double(err))>tol
        error("case_2D_twisted failed: Part 1 (scalar, single-group).")
    end

    %% Part 2: single-group, VECTOR state (dims=2), degree 2, twisted order
    pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
    Kvec = [s1_dum1*s2_dum1^2, 5*s1_dum2*s2_dum2^3, -2*s1_dum2^2*s2_dum1, s1_dum1+s2_dum2];
    KopB = intop({Kvec},intmat,vars,dom,false);
    KopB2 = combine_terms(KopB,2,2,depmat);

    vxvec = [s^2*t+1; 3*t^3-s*t^2+2];
    v1vec = subs(vxvec,[s;t],[s1_1;s2_1]);
    v2vec = subs(vxvec,[s;t],[s1_2;s2_2]);
    xvalB = kron(v1vec,v2vec);

    Y0B = apply_op(KopB,xvalB);
    Y1B = apply_op(KopB2,xvalB);
    errB = cleanpoly(Y0B-Y1B,1e-10);
    if ~isdouble(errB) || abs(double(errB))>tol
        error("case_2D_twisted failed: Part 2 (vector, single-group).")
    end

    %% Part 3: single-group, SCALAR state, degree 3, twisted order
    pvar s1_dum1 s1_dum2 s1_dum3 s2_dum1 s2_dum2 s2_dum3
    I1c = [0,3, 0,1, 0,4];   % s1: slot2 inner[0,slot1], slot1 mid[0,slot3], slot3 outer
    I2c = [0,4, 0,1, 0,2];   % s2: slot1 outer, slot2 mid[0,slot1], slot3 inner[0,slot2]
    intmatC = {I1c,I2c};
    Kc = s1_dum1*s2_dum1 + 2*s1_dum2*s2_dum2^2 - 3*s1_dum3^2*s2_dum3 + s1_dum1*s1_dum2*s2_dum3;
    KopC = intop({Kc},intmatC,vars,dom,false);
    KopC2 = combine_terms(KopC,3,1,depmat);

    pvar s1_3 s2_3
    vxC1 = subs(vx,[s;t],[s1_1;s2_1]);
    vxC2 = subs(vx,[s;t],[s1_2;s2_2]);
    vxC3 = subs(vx,[s;t],[s1_3;s2_3]);
    xvalC = vxC1*vxC2*vxC3;

    Y0C = apply_op(KopC,xvalC);
    Y1C = apply_op(KopC2,xvalC);
    errC = cleanpoly(Y0C-Y1C,1e-10);
    if ~isdouble(errC) || abs(double(errC))>tol
        error("case_2D_twisted failed: Part 3 (scalar, degree 3, single-group).")
    end

    %% Part 4: MULTI-group, SCALAR state -- the key stress test. Two terms
    % along s1 with OPPOSITE nesting order of the same factor-slot pair
    % (forcing combine_terms to split into 2 groups, apply a DIFFERENT
    % canonicalizing permutation to each, then sum), plus a twisted order
    % along s2.
    I1m = [0,2,0,3;  0,3,0,1];   % s1: term1 slot1-inner/slot2-outer, term2 REVERSED
    I2m = [0,3,0,1];             % s2: slot1 outer, slot2 inner (one term)
    intmatM = {I1m,I2m};
    pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
    K1m = s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*s1_dum2^2*s2_dum1;
    K2m = 3*s1_dum1^2*s2_dum2 - s1_dum2*s2_dum1^2 + 7;
    KopM = intop({K1m;K2m},intmatM,vars,dom,false);
    KopM2 = combine_terms(KopM,2,1,depmat);
    % Confirm this genuinely exercised the multi-group path (not a no-op):
    if size(KopM.intmat{1},1)<2
        error("case_2D_twisted failed: Part 4 setup assumption violated (need 2+ terms along s1).")
    end

    xvalM = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);
    Y0M = apply_op(KopM,xvalM);
    Y1M = apply_op(KopM2,xvalM);
    errM = cleanpoly(Y0M-Y1M,1e-10);
    if ~isdouble(errM) || abs(double(errM))>tol
        error("case_2D_twisted failed: Part 4 (scalar, MULTI-group, twisted).")
    end

    %% Part 5: MULTI-group, VECTOR state (dims=2), same twisted multi-group setup.
    pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
    K1mv = [s1_dum1*s2_dum1^2, 5*s1_dum2*s2_dum2^3, -2*s1_dum2^2*s2_dum1, s1_dum1+s2_dum2];
    K2mv = [3*s1_dum1^2*s2_dum2, -s1_dum2*s2_dum1^2+7, s1_dum1-s1_dum2, s2_dum1*s2_dum2];
    KopMv = intop({K1mv;K2mv},intmatM,vars,dom,false);
    KopMv2 = combine_terms(KopMv,2,2,depmat);

    xvalMv = kron(subs(vxvec,[s;t],[s1_1;s2_1]), subs(vxvec,[s;t],[s1_2;s2_2]));
    Y0Mv = apply_op(KopMv,xvalMv);
    Y1Mv = apply_op(KopMv2,xvalMv);
    errMv = cleanpoly(Y0Mv-Y1Mv,1e-10);
    if ~isdouble(errMv) || abs(double(errMv))>tol
        error("case_2D_twisted failed: Part 5 (vector, MULTI-group, twisted).")
    end

    disp('case_2D_twisted passed.')
end

%% --------------------------------------------------------------------------
function case_canonical()
% case_canonical -- Tests that @intop/combine_terms.m
% (the form with arguments (Kop,degs,dims,depmat)) gives the same intmat
% and the same kernels for equal functionals, so that the coefficients of
% equal functionals can be compared cell by cell (as piesos_eq does).
%
% Every integral is written as a sum of cells, one for each ordering of the
% dummy variables consistent with its limits, in canonical form (each
% variable integrated from a up to the next larger one, the largest over
% the full domain, multipliers [t,t] at a variable that is not a
% multiplier). For copies of the same state, the cells are then relabeled
% to the same canonical order, adding the kernels of cells that coincide.
% This is checked
%   Part 1: for the square [0,3,0,3] with kernel k, which is the chain
%           [0,3,0,1] with k(s1,s2)+k(s2,s1), and for other forms of the
%           same region;
%   Part 2: for two different states, for which the cells cannot be
%           relabeled;
%   Part 3: for the cube [0,4,0,4,0,4], which becomes one chain, with the
%           kernel the sum over the 6 orderings;
%   Part 4: for multipliers, whose kernels are written in terms of the
%           variable at which they are evaluated: the kernel of [0,3,1,1]
%           is k(s1,s1);
%   Part 5: for random functionals with ties in N=1,2,3 spatial variables:
%           relabeling the copies of the functional (reorder_states) does
%           not change the result of combine_terms, the result is made of
%           cells, the value is unchanged, and calling combine_terms again
%           does nothing;
%   Part 6: the same, with multipliers, for kernels that depend on all
%           variables, also those of the multipliers, which are replaced
%           by the variable at which the multiplier is evaluated (the
%           lowest of the equal variables);
%   Part 7: for rows that are not ordered, e.g. s1 integrated from s3 to s2,
%           which the 'intop' constructor rejects: after splitting them with
%           canonicalize_intmat (the indicator of s3<=s1<=s2; an empty region
%           gives an error), the value is unchanged and the result is made of
%           cells;
%   Part 8: for groups of equal variables (multipliers), also for
%           different states;
%   Part 9: for vector-valued states: copies that are equal in all spatial
%           variables are interchangeable, x_i(s)*x_j(s) = x_j(s)*x_i(s),
%           so that only the sum of the kernels on the columns (i,j) and
%           (j,i) is determined, which is put on the column (i,j), i<=j;
%   Part 10: for kernels with decision variables, as the unknown
%           functionals of the SOS programs.
%

% NOTES
% - Oracle: expected intmat and kernels are written out by hand per case;
%   the value check compares apply_op before vs. after combine_terms.
% - Random parts (5, 6, 9b) are seeded (rng(60000+..), rng(61000+..)).

    tol = 1e-8;
    pvar s_dum1 s_dum2 s_dum3
    sw = @(k,a,b) subs(k,[a;b],[b;a]);

    %% Part 1: d=2, one variable: the square and its chains
    k = s_dum1^2+3*s_dum1*s_dum2^3;                 % not symmetric
    ksw = sw(k,s_dum1,s_dum2);
    mk = @(K,im) intop({K},{im},{'s'},[0,1],false);
    % intmat, kernel, expected kernel on [0,3,0,1]
    cases = {[0,3,0,3],k,k+ksw; [0,3,0,3],ksw,k+ksw; [0,3,0,1],k+ksw,k+ksw; [0,2,0,3],k+ksw,k+ksw; ...
             [2,3,0,3],k,k; [0,3,0,1],k,k; [0,2,0,3],k,ksw};
    for j=1:size(cases,1)
        K1 = combine_terms(mk(cases{j,2},cases{j,1}),2,1,true);
        if ~isequal(K1.intmat{1},[0,3,0,1])
            error("case_canonical failed: Part 1, case %d: intmat %s, expected [0,3,0,1].",j,mat2str(K1.intmat{1}))
        end
        kk = kernels(K1);
        err = cleanpoly(kk{1}-cases{j,3},1e-12);
        if ~isdouble(err) || abs(double(err))>tol
            error("case_canonical failed: Part 1, case %d: wrong kernel.",j)
        end
        if ~same_op(combine_terms(K1,2,1,true),K1)
            error("case_canonical failed: Part 1, case %d: not idempotent.",j)
        end
    end

    %% Part 2: two different states: [0,3,0,3] is the sum of two cells, with the same kernel
    K1 = combine_terms(mk(k,[0,3,0,3]),[1,1],[1,1],true(2,1));
    kk = kernels(K1);
    if ~isequal(K1.intmat{1},[0,2,0,3;0,3,0,1]) && ~isequal(K1.intmat{1},[0,3,0,1;0,2,0,3])
        error("case_canonical failed: Part 2, intmat %s.",mat2str(K1.intmat{1}))
    end
    for r=1:2
        err = cleanpoly(kk{r}-k,1e-12);
        if ~isdouble(err) || abs(double(err))>tol
            error("case_canonical failed: Part 2, kernel %d.",r)
        end
    end
    K1 = combine_terms(mk(k,[2,3,0,3]),[1,1],[1,1],true(2,1));
    if ~isequal(K1.intmat{1},[0,3,0,1])
        error("case_canonical failed: Part 2, [2,3,0,3] is not [0,3,0,1].")
    end

    %% Part 3: the cube is one chain, with the kernel summed over all orderings
    k3 = s_dum1+2*s_dum2^2+3*s_dum3^3;
    K1 = combine_terms(intop({k3},{[0,4,0,4,0,4]},{'s'},[0,1],false),3,1,true);
    if ~isequal(K1.intmat{1},[0,4,0,1,0,2])
        error("case_canonical failed: Part 3, intmat %s.",mat2str(K1.intmat{1}))
    end
    ksum = 0;
    P = perms(1:3);
    vv = [s_dum1;s_dum2;s_dum3];
    for r=1:6
        ksum = ksum+subs(k3,vv,vv(P(r,:)));
    end
    kk = kernels(K1);
    err = cleanpoly(kk{1}-ksum,1e-12);
    if ~isdouble(err) || abs(double(err))>tol
        error("case_canonical failed: Part 3, kernel is not the sum over the orderings.")
    end

    %% Part 4: multipliers: [2,2,0,3] is [0,3,1,1] with the variables switched
    K1 = combine_terms(mk(k,[0,3,1,1]),2,1,true);
    K2 = combine_terms(mk(k,[2,2,0,3]),2,1,true);
    if ~isequal(K1.intmat{1},[0,3,1,1]) || ~isequal(K2.intmat{1},[0,3,1,1])
        error("case_canonical failed: Part 4, intmat %s and %s, expected [0,3,1,1].",mat2str(K1.intmat{1}),mat2str(K2.intmat{1}))
    end
    % Both are evaluated at s2=s1: the kernel is k(s1,s1), without s_dum2
    kd = subs(k,s_dum2,s_dum1);
    kk1 = kernels(K1);      kk2 = kernels(K2);
    err1 = cleanpoly(kk1{1}-kd,1e-12);     err2 = cleanpoly(kk2{1}-kd,1e-12);
    if ~isdouble(err1) || abs(double(err1))>tol || ~isdouble(err2) || abs(double(err2))>tol
        error("case_canonical failed: Part 4, kernels.")
    end

    %% Part 5 and 6: random functionals, N in {1,2,3}, with ties; without and with multipliers
    for variant = 1:2
        pdelta = 0.3*(variant==2);
        for trial = 1:25
            rng(60000+100*variant+trial);
            N = randi(3);
            d = randi([2,4]);
            if N>=2,    d = min(d,3);   end
            nterms_i = randi(3,1,N);
            vnames = arrayfun(@(a) ['s',num2str(a)],1:N,'UniformOutput',false);
            vars = struct('in',{vnames},'out',{cell(1,0)});
            dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
            intmat = cell(1,N);
            dumvars = cell(1,N);
            for i=1:N
                intmat{i} = rand_ties_intmat(d,nterms_i(i),pdelta);
                dumvars{i} = polynomial(arrayfun(@(a) [vnames{i},'_dum',num2str(a)],(1:d)','UniformOutput',false));
            end
            % Kernels in all dummy variables, also those of multipliers
            allnames = cell(1,0);
            for i=1:N
                allnames = [allnames,dumvars{i}.varname(:)'];   %#ok<AGROW>
            end
            params = cell([nterms_i,1]);
            for lin=1:numel(params)
                params{lin} = rand_poly([1,1],allnames,1,6);
            end
            Kop = intop(params,intmat,vars,dom,false);
            K1 = combine_terms(Kop,d,1,true(1,N));
            % (a) idempotent
            if ~same_op(combine_terms(K1,d,1,true(1,N)),K1)
                error("case_canonical failed: Part %d, trial %d (N=%d,d=%d): not idempotent.",4+variant,trial,N,d)
            end
            % (b) all integrals are cells
            for i=1:N
                for r=1:size(K1.intmat{i},1)
                    if ~is_cell_row(K1.intmat{i}(r,:))
                        error("case_canonical failed: Part %d, trial %d: row %s is not a cell.",4+variant,trial,mat2str(K1.intmat{i}(r,:)))
                    end
                end
            end
            % (b2) the kernels do not depend on the variables of multipliers
            kk = kernels(K1);
            for lin=1:numel(kk)
                kc = cleanpoly(kk{lin},1e-12);      % a zero kernel keeps the names of its variables
                sb = cell(1,N);
                [sb{:}] = ind2sub([size(K1.params,1:N),1],lin);
                for i=1:N
                    row = K1.intmat{i}(sb{i},:);
                    is_dlt = row(1:2:end)==row(2:2:end);
                    if any(ismember(kc.varname,dumvars{i}.varname(is_dlt(:))))
                        error("case_canonical failed: Part %d, trial %d (N=%d,d=%d): a kernel depends on a multiplier variable.",4+variant,trial,N,d)
                    end
                end
            end
            % (c) the value is unchanged
            tvars = polynomial(arrayfun(@(a) ['t',num2str(a)],(1:N)','UniformOutput',false));
            vx = rand_poly([1,1],tvars,2,6);
            if all(isequal(vx,0)),  vx = polynomial(1);     end
            xval = 1;
            for kk=1:d
                copyvars = arrayfun(@(a) [vnames{a},'_',num2str(kk)],(1:N)','UniformOutput',false);
                xval = kron(xval,subs(vx,tvars,polynomial(copyvars)));
            end
            Y0 = apply_op(Kop,xval);
            err = cleanpoly(Y0-apply_op(K1,xval),1e-10);
            if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,norm(double(cleanpoly(Y0,1e-10)),'fro'))
                error("case_canonical failed: Part %d, trial %d (N=%d,d=%d): the functional changed.",4+variant,trial,N,d)
            end
            % (d) relabeling the copies does not change the result
            perm = randperm(d);
            K2 = combine_terms(reorder_states(Kop,perm,ones(1,d),true(d,N)),d,1,true(1,N));
            if ~same_kernels(K1,K2)
                error("case_canonical failed: Part %d, trial %d (N=%d,d=%d): relabeled copies give another form.",4+variant,trial,N,d)
            end
        end
    end

    %% Part 7: rows that are not ordered: [3,2,0,4,0,4] is s1 from s3 to s2, nonzero only for s3<=s1<=s2
    pvar s s_1 s_2 s_3
    vb = {1+s+2*s^2, 2-s+s^3};
    rows7 = {[3,2,0,3,0,4],[2,3,0,4,0,4],[3,2,0,4,0,4],[0,4,3,1,0,4]};
    degs7 = {3,3,3,[2,1]};
    for j=1:numel(rows7)
        for r=1:numel(degs7)
            dg = degs7{r};
            if numel(dg)==1,    state_idcs = [1,1,1];   else,   state_idcs = [1,1,2];  end
            k7 = rand_poly([1,1],{'s_dum1','s_dum2','s_dum3'},1,5);
            % The constructor rejects rows that are not ordered
            try
                [im7,p7] = intop.canonicalize_intmat(rows7(j),{k7});
            catch
                % [3,2,0,3,0,4]: s2<=s3, so that s3<=s1<=s2 only on a null set: the functional is zero
                if j~=1,    error('case_canonical failed: Part 7, row %s: empty region.',mat2str(rows7{j})),   end
                continue
            end
            Kop = intop(p7,im7,{'s'},[0,1],false);
            K1 = combine_terms(Kop,dg,ones(1,numel(dg)),true(numel(dg),1));
            if ~same_op(combine_terms(K1,dg,ones(1,numel(dg)),true(numel(dg),1)),K1)
                error("case_canonical failed: Part 7, row %s, degs %s: not idempotent.",mat2str(rows7{j}),mat2str(dg))
            end
            for q=1:size(K1.intmat{1},1)
                if ~is_cell_row(K1.intmat{1}(q,:))
                    error("case_canonical failed: Part 7, row %s: %s is not a cell.",mat2str(rows7{j}),mat2str(K1.intmat{1}(q,:)))
                end
            end
            sv = [s_1,s_2,s_3];
            xval = 1;
            for c=1:3
                xval = kron(xval,subs(vb{state_idcs(c)},s,sv(c)));
            end
            Y0 = apply_op(Kop,xval);
            err = cleanpoly(Y0-apply_op(K1,xval),1e-10);
            if ~isdouble(err) || abs(double(err))>tol*max(1,abs(double(cleanpoly(Y0,1e-10))))
                error("case_canonical failed: Part 7, row %s, degs %s: the functional changed.",mat2str(rows7{j}),mat2str(dg))
            end
        end
    end
    %% Part 8: groups of equal variables: the lowest variable of the group is free, the others are evaluated at it
    % [0,4,1,1,1,1]: s2=s3=s1; [3,3,3,3,0,4]: s1=s3, s2=s3, s3 free (the same group, other representative)
    k3 = s_dum1+2*s_dum2^2+3*s_dum3^3;
    kd3 = subs(k3,[s_dum2;s_dum3],[s_dum1;s_dum1]);             % s1+2*s1^2+3*s1^3
    forms8 = {[0,4,1,1,1,1],[3,3,3,3,0,4],[0,4,1,1,2,2],[0,4,3,3,1,1]};
    for j=1:numel(forms8)
        for r=1:2
            if r==1,    dg = 3;         dm = 1;         dp = true;
            else,       dg = [1,1,1];   dm = [1,1,1];   dp = true(3,1);    end
            K1 = combine_terms(intop({k3},{forms8{j}},{'s'},[0,1],false),dg,dm,dp);
            if ~isequal(K1.intmat{1},[0,4,1,1,1,1])
                error("case_canonical failed: Part 8, %s with degs %s: intmat %s, expected [0,4,1,1,1,1].",mat2str(forms8{j}),mat2str(dg),mat2str(K1.intmat{1}))
            end
            kk = kernels(K1);
            err = cleanpoly(kk{1}-kd3,1e-12);
            if ~isdouble(err) || abs(double(err))>tol
                error("case_canonical failed: Part 8, %s with degs %s: kernel is not k(s1,s1,s1).",mat2str(forms8{j}),mat2str(dg))
            end
            if ~same_op(combine_terms(K1,dg,dm,dp),K1)
                error("case_canonical failed: Part 8, %s: not idempotent.",mat2str(forms8{j}))
            end
        end
    end
    % two groups, d=4: [2,2,0,4,4,4,0,5] is s1=s2 (below the group s3=s4, s4 free): free s1 and s3
    pvar s_dum4
    k4 = s_dum1+2*s_dum2^2+3*s_dum3^3+s_dum4*s_dum2;
    kd4 = subs(k4,[s_dum2;s_dum4],[s_dum1;s_dum3]);
    K1 = combine_terms(intop({k4},{[2,2,0,4,4,4,0,5]},{'s'},[0,1],false),[1,1,1,1],[1,1,1,1],true(4,1));
    kk = kernels(K1);
    err = cleanpoly(kk{1}-kd4,1e-12);
    if ~isequal(K1.intmat{1},[0,3,1,1,0,5,3,3]) || ~isdouble(err) || abs(double(err))>tol
        error("case_canonical failed: Part 8, two groups: intmat %s.",mat2str(K1.intmat{1}))
    end
    %% Part 9: vector-valued state (2 elements), copies equal at a point: the columns (1,2) and (2,1) are added
    % kernel columns in the order (1,1),(1,2),(2,1),(2,2)
    kc = {s_dum1+1, 2*s_dum1^2, -3*s_dum1, 5};
    K9 = intop({[kc{1},kc{2},kc{3},kc{4}]},{[0,3,1,1]},{'s'},[0,1],false);
    K1 = combine_terms(K9,2,2,true);
    kk = kernels(K1);
    kexp = [kc{1},kc{2}+kc{3},0,kc{4}];
    err = cleanpoly(kk{1}-kexp,1e-12);
    if ~isequal(K1.intmat{1},[0,3,1,1]) || ~isdouble(err) || any(abs(double(err))>tol)
        error("case_canonical failed: Part 9, the columns of equal copies are not combined.")
    end
    if ~same_op(combine_terms(K1,2,2,true),K1)
        error("case_canonical failed: Part 9, not idempotent.")
    end
    % the value does not change, and the symmetric and asymmetric kernels of the same functional agree
    pvar s s_1 s_2
    vx = [1+s; 2-s^2];
    xval = kron(subs(vx,s,s_1),subs(vx,s,s_2));
    % evaluated at s_2=s_1 through the multiplier: apply_op does it
    Y0 = apply_op(K9,xval);     Y1 = apply_op(K1,xval);
    kc_sym = {kc{1},kc{2}+kc{3},0,kc{4}};
    K9s = intop({[kc_sym{1},kc_sym{2},kc_sym{3},kc_sym{4}]},{[0,3,1,1]},{'s'},[0,1],false);
    Y2 = apply_op(K9s,xval);
    if ~isdouble(cleanpoly(Y0-Y1,1e-10)) || ~isdouble(cleanpoly(Y0-Y2,1e-10)) || any(abs(double(cleanpoly(Y0-Y1,1e-10)))>tol) || any(abs(double(cleanpoly(Y0-Y2,1e-10)))>tol)
        error("case_canonical failed: Part 9, the functional changed.")
    end
    K2 = combine_terms(K9s,2,2,true);
    if ~same_op(K1,K2)
        error("case_canonical failed: Part 9, equal functionals have different kernels.")
    end

    %% Part 9b: random vector-valued functionals with multipliers: idempotent, the value is unchanged, relabeling the copies does not change the result
    for trial = 1:15
        rng(61000+trial);
        N = randi(2);
        d = randi([2,3]);
        nterms_i = randi(3,1,N);
        vnames = arrayfun(@(a) ['s',num2str(a)],1:N,'UniformOutput',false);
        vars = struct('in',{vnames},'out',{cell(1,0)});
        dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
        intmat = cell(1,N);
        allnames = cell(1,0);
        for i=1:N
            intmat{i} = rand_ties_intmat(d,nterms_i(i),0.4);
            allnames = [allnames,arrayfun(@(a) [vnames{i},'_dum',num2str(a)],1:d,'UniformOutput',false)];  %#ok<AGROW>
        end
        params = cell([nterms_i,1]);
        for lin=1:numel(params)
            params{lin} = rand_poly([1,2^d],allnames,1,5);
        end
        Kop = intop(params,intmat,vars,dom,false);
        K1 = combine_terms(Kop,d,2,true(1,N));
        if ~same_op(combine_terms(K1,d,2,true(1,N)),K1)
            error("case_canonical failed: Part 9b, trial %d (N=%d,d=%d): not idempotent.",trial,N,d)
        end
        tvars = polynomial(arrayfun(@(a) ['t',num2str(a)],(1:N)','UniformOutput',false));
        vx = rand_poly([2,1],tvars,2,5);
        xval = 1;
        for kcopy=1:d
            copyvars = arrayfun(@(a) [vnames{a},'_',num2str(kcopy)],(1:N)','UniformOutput',false);
            xval = kron(xval,subs(vx,tvars,polynomial(copyvars)));
        end
        Y0 = apply_op(Kop,xval);
        err = cleanpoly(Y0-apply_op(K1,xval),1e-10);
        if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,norm(double(cleanpoly(Y0,1e-10)),'fro'))
            error("case_canonical failed: Part 9b, trial %d (N=%d,d=%d): the functional changed.",trial,N,d)
        end
        perm = randperm(d);
        K2 = combine_terms(reorder_states(Kop,perm,2*ones(1,d),true(d,N)),d,2,true(1,N));
        if ~same_kernels(K1,K2)
            error("case_canonical failed: Part 9b, trial %d (N=%d,d=%d): relabeled copies give another form.",trial,N,d)
        end
    end
    %% Part 10: decision variables in the kernel of a multiplier
    dpvar c1 c2
    Kd = intop({c1*s_dum1+c2*s_dum2^2+3},{[0,3,1,1]},{'s'},[0,1],false);
    Kc = combine_terms(Kd,2,1,true);
    kk = kernels(Kc);
    ed = kk{1}-(c1*s_dum1+c2*s_dum1^2+3);
    if ~isequal(Kc.intmat{1},[0,3,1,1]) || (~isempty(ed.C) && max(abs(ed.C(:)))>tol)
        error("case_canonical failed: Part 10, kernel with decision variables on a multiplier.")
    end
    if ~same_op(combine_terms(Kc,2,1,true),Kc)
        error("case_canonical failed: Part 10, not idempotent.")
    end
    disp('case_canonical passed.')
end

%% --------------------------------------------------------------------------
function case_idempotent()
% case_idempotent -- Tests that
% @intop/combine_terms.m (the form with arguments (Kop,degs,dims,depmat))
% puts the repeated factors of a state in a canonical order, so that
%   (1) combining the terms of a functional that was already combined
%       leaves it unchanged, exactly: the same intmat, monomials, and
%       parameters;
%   (2) the combined functional represents the same functional: it gives
%       the same value on x o ... o x as the original;
%   (3) in every term with nonzero parameters, the factors are ordered by
%       decreasing order of integration (the last variable integrated
%       first), ties broken by the order of integration along the next
%       spatial variable, and remaining ties keeping their order. This is
%       checked independently of the code, from the integration orders
%       (int_order) of the combined functional. Terms without parameters
%       are skipped: the terms are a product over the rows along each
%       variable, so they cannot all be in order, and combine_terms leaves
%       them as they are.
%
%
% NOTES
% - Part 1 expects the tied square [0,3,0,3] as the chain [0,3,0,1] with
%   kernel k(s1,s2)+k(s2,s1): combine_terms splits every integral into
%   cells (see Test_intop_combine_terms/case_canonical).
% - Part 3 draws tied orders of integration with rand_ties_intmat.m (seeds
%   50000+trial).

    tol = 1e-8;

    %% Part 1: N=1, degree 2: tied full square, and chains (ordered)
    pvar s_dum1 s_dum2
    intmat_arr = {[0,3,0,3],[0,3,0,1],[0,2,0,3]};
    for j=1:numel(intmat_arr)
        Kop = intop({s_dum1*(1+s_dum2^2)},intmat_arr(j),{'s'},[0,1],false);
        K1 = combine_terms(Kop,2,1,true);
        K2 = combine_terms(K1,2,1,true);
        K3 = combine_terms(K2,2,1,true);
        if ~same_op_idempotent(K1,K2) || ~same_op_idempotent(K2,K3)
            error("case_idempotent failed: Part 1, intmat %s is not left unchanged by a second call.",mat2str(intmat_arr{j}))
        end
        % Expected kernel: k for [0,3,0,1], k(s2,s1) for [0,2,0,3], k(s1,s2)+k(s2,s1) for the square
        if ~isequal(K1.intmat{1},[0,3,0,1])
            error("case_idempotent failed: Part 1, intmat %s not given as [0,3,0,1].",mat2str(intmat_arr{j}))
        end
        ksw = s_dum2*(1+s_dum1^2);
        kexp = {s_dum1*(1+s_dum2^2)+ksw, s_dum1*(1+s_dum2^2), ksw};
        kk = kernels(K1);
        err = cleanpoly(kk{1}-kexp{j},1e-12);
        if ~isdouble(err) || abs(double(err))>tol
            error("case_idempotent failed: Part 1, kernel of intmat %s.",mat2str(intmat_arr{j}))
        end
    end

    %% Part 2: N=2, degree 2, hand-picked cases
    vars = struct('in',{{'s1','s2'}},'out',{cell(1,0)});
    dom = struct('in',[0,1;0,1],'out',zeros(0,2));
    depmat = true(1,2);
    pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
    K1m = s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*s1_dum2^2*s2_dum1;
    K2m = 3*s1_dum1^2*s2_dum2 - s1_dum2*s2_dum1^2 + 7;
    K3m = s1_dum1 + 2*s2_dum2 + s1_dum2*s2_dum1;
    % intmat along s1, along s2, kernels (one per combination of terms)
    cases = {{[0,2,0,3;0,3,0,1],[0,3,0,1],{K1m;K2m}}, ...              % 2 orders along s1, twisted along s2
             {[0,2,0,3;0,3,0,1;0,3,0,3],[0,3,0,1],{K1m;K2m;K3m}}, ...  % plus a tied term along s1
             {[0,3,0,3],[0,3,0,1;0,2,0,3],{K1m,K2m}}, ...              % tied along s1, orders along s2
             {[0,3,0,3],[0,3,0,3],{K3m}}, ...                          % tied along both
             {[0,2,0,3],[0,2,0,3;0,3,0,1],{K1m,K2m}}};                 % ordered along both, opposite along s2

    pvar s t s1_1 s2_1 s1_2 s2_2
    vx = s^2*t + 3*t^3 - s*t^2 + 1;
    xval = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);
    for j=1:numel(cases)
        Kop = intop(cases{j}{3},cases{j}(1:2),vars,dom,false);
        K1 = combine_terms(Kop,2,1,depmat);
        K2 = combine_terms(K1,2,1,depmat);
        K3 = combine_terms(K2,2,1,depmat);
        if ~same_op_idempotent(K1,K2) || ~same_op_idempotent(K2,K3)
            error("case_idempotent failed: Part 2, case %d is not left unchanged by a second call.",j)
        end
        if ~is_canon(K1,2)
            error("case_idempotent failed: Part 2, case %d: factors not in canonical order.",j)
        end
        err = cleanpoly(apply_op(Kop,xval)-apply_op(K1,xval),1e-10);
        if ~isdouble(err) || abs(double(err))>tol
            error("case_idempotent failed: Part 2, case %d: the functional changed.",j)
        end
    end

    %% Part 3: random N in {1,2,3}, degree in {2,3,4}, with ties, scalar or vector state
    for trial = 1:40
        rng(50000+trial);
        N = randi(3);
        d = randi([2,4]);
        if N==3,    d = min(d,3);   end
        if N==2,    d = min(d,3);   end
        dimval = randi(2);
        if d>=3,    dimval = 1;     end
        nterms_i = randi(3,1,N);

        vnames = cell(1,N);
        for i=1:N,  vnames{i} = ['s',num2str(i)];   end
        vars = struct('in',{vnames},'out',{cell(1,0)});
        dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
        intmat = cell(1,N);
        dumvars = cell(1,N);
        for i=1:N
            intmat{i} = rand_ties_intmat(d,nterms_i(i),0.15);
            dv = cell(d,1);
            for j=1:d,  dv{j} = [vnames{i},'_dum',num2str(j)];  end
            dumvars{i} = polynomial(dv);
        end
        allpvars = vertcat(dumvars{:});
        params = cell([nterms_i,1]);
        for lin=1:numel(params)
            params{lin} = rand_poly([1,dimval^d],allpvars,1,6);
        end
        Kop = intop(params,intmat,vars,dom,false);
        K1 = combine_terms(Kop,d,dimval,true(1,N));
        K2 = combine_terms(K1,d,dimval,true(1,N));
        if ~same_op_idempotent(K1,K2)
            error("case_idempotent failed: Part 3, trial=%d, N=%d, d=%d: not left unchanged by a second call.",trial,N,d)
        end
        if ~is_canon(K1,d)
            error("case_idempotent failed: Part 3, trial=%d, N=%d, d=%d: factors not in canonical order.",trial,N,d)
        end

        % Same-state input x^{od}: one base formula, renamed per copy
        tvars_p = polynomial(cellfun(@(a) ['t',num2str(a)],num2cell(1:N)','UniformOutput',false));
        vx = rand_poly([dimval,1],tvars_p,2,6);
        if all(isequal(vx,0)),  vx = ones(dimval,1);    end
        xval = 1;
        for k=1:d
            copyvars = cellfun(@(a) [a,'_',num2str(k)],vnames(:),'UniformOutput',false);
            xval = kron(xval,subs(vx,tvars_p,polynomial(copyvars)));
        end
        Y0 = apply_op(Kop,xval);
        Y1 = apply_op(K1,xval);
        err = cleanpoly(Y0-Y1,1e-10);
        nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
        if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
            error("case_idempotent failed: Part 3, trial=%d, N=%d, d=%d: the functional changed.",trial,N,d)
        end
    end

    %% Part 4: decision variables in the kernels (dpvar), degree 3, one variable, tied and ordered terms
    pvar s_dum1 s_dum2 s_dum3
    dpvar c1 c2
    Kd = {c1*s_dum1+2*s_dum3; c2*s_dum2^2+s_dum1*s_dum3; c1*s_dum3+c2};
    Kop = intop(Kd,{[0,4,0,4,0,4;0,2,0,3,0,4;0,4,0,1,0,2]},{'s'},[0,1],false);
    K1 = combine_terms(Kop,3,1,true);
    K2 = combine_terms(K1,3,1,true);
    if ~same_op_idempotent(K1,K2) || ~is_canon(K1,3)
        error("case_idempotent failed: Part 4, kernels with decision variables.")
    end
    pvar s s_1 s_2 s_3
    vb = 1+s+2*s^2;
    xval = subs(vb,s,s_1)*subs(vb,s,s_2)*subs(vb,s,s_3);
    errd = apply_op(Kop,xval)-apply_op(K1,xval);
    if ~isempty(errd.C) && max(abs(errd.C(:)))>tol
        error("case_idempotent failed: Part 4, the functional with decision variables changed.")
    end

    disp('case_idempotent passed.')
end

%% --------------------------------------------------------------------------
function case_ND_random()
% case_ND_random -- Randomized test of
% @intop/combine_terms.m for a repeated state spanning N>=2 spatial
% directions, complementing the hand-picked scenarios in
% Test_intop_combine_terms/case_2D_twisted.
%
% Each trial draws N in {2,3} spatial directions, degree d in {2,3} (number
% of copies of the single state), state size dims in {1,2}, and 1-3 terms
% per direction, each term's nested integration order drawn independently
% at random (rand_valid_intmat.m). Verifies that
% apply_op(Kop,xval) == apply_op(combine_terms(Kop,...),xval) for a
% same-state input (one base formula, renamed per copy), which is the input
% family combine_terms is used on.
%
% Oracle: combine_terms only permutes copies of the same state
%   (reorder_states.m, Kop_out*(z_{so(1)} o ... o z_{so(d)}) =
%   Kop*(z_1 o ... o z_d) for any z_i and permutation 'so', see
%   Test_intop_reorder_states.m). With z_1=...=z_d=x the input is invariant
%   under permutation of its factors, so the value on x^{od} is unchanged
%   however the terms are grouped or reordered.
%
%
% NOTES
% - Value check only: it does not check the canonical form (see
%   Test_intop_combine_terms/case_canonical).

    tol = 1e-8;

    for trial = 1:20
        rng(30000+trial);
        N = randi([2,3]);
        d = randi([2,3]);
        dimval = randi([1,2]);
        if d==3 && N==3
            dimval = 1;   % keep input dimension n=dims^d manageable
        end
        nterms_i = randi([1,3],1,N);

        vnames = cell(1,N);
        for i=1:N, vnames{i} = ['s',num2str(i)]; end
        vars = struct('in',{vnames},'out',{cell(1,0)});
        dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));

        intmat = cell(1,N);
        dumvars = cell(1,N);
        for i=1:N
            intmat{i} = rand_valid_intmat(d,nterms_i(i));
            dv = cell(d,1);
            for j=1:d, dv{j} = [vnames{i},'_dum',num2str(j)]; end
            dumvars{i} = polynomial(dv);
        end
        allpvars = vertcat(dumvars{:});

        n_in = dimval^d;
        params = cell(nterms_i);
        for lin=1:numel(params)
            params{lin} = rand_poly([1,n_in],allpvars,2,8);
        end

        Kop = intop(params,intmat,vars,dom,false);
        Kop2 = combine_terms(Kop,d,dimval,true(1,N));

        % Same-state input x^{od}: one base formula, renamed per copy
        tvars = cell(1,N);
        for i=1:N, tvars{i} = ['t',num2str(i)]; end
        tvars_p = polynomial(tvars(:));
        vx = rand_poly([dimval,1],tvars_p,2,6);
        if all(isequal(vx,0)), vx = ones(dimval,1); end

        xval = 1;
        for k=1:d
            copyvars = cell(1,N);
            for i=1:N, copyvars{i} = [vnames{i},'_',num2str(k)]; end
            vcopy = subs(vx,tvars_p,polynomial(copyvars(:)));
            xval = kron(xval,vcopy);
        end

        Y0 = apply_op(Kop,xval);
        Y1 = apply_op(Kop2,xval);
        while isa(Y0,'cell'), Y0 = Y0{1}; end
        while isa(Y1,'cell'), Y1 = Y1{1}; end
        err = cleanpoly(Y0-Y1,1e-10);
        nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
        if ~isdouble(err) || norm(double(err),'fro') > tol*max(1,nrmref)
            error("case_ND_random failed: trial=%d, N=%d, d=%d, dims=%d, nterms=%s.", ...
                trial,N,d,dimval,mat2str(nterms_i))
        end
    end

    disp('case_ND_random passed.')
end

%% --------------------------------------------------------------------------
function case_sum_blocks()
% case_sum_blocks -- Tests that
% @intop/combine_terms.m, which sums the blocks of terms with the same
% order of the copies of a state in one pass, gives a functional that
%   (1) represents the same functional as the input: the same value on
%       x o ... o x, for kernels that are polynomials or depend on decision
%       variables (dpvar);
%   (2) keeps the decision variables Zd of the input (the blocks share
%       them, they are not compared or merged);
%   (3) has unique, sorted integrals, and is left unchanged by a second
%       call;
%   (4) does not depend on the order of the terms of the input: permuting
%       the rows of intmat (and the parameters accordingly) gives the same
%       functional (the same value, and the same intmat after combining).
%
%
% NOTES
% - Does not compare the one-pass sum of the blocks with a sum of the
%   blocks by plus; it checks only the properties (1)-(4).

    tol = 1e-8;

    %% Part 1: N=2, degree 2, twisted terms with decision variables
    vars = struct('in',{{'s1','s2'}},'out',{cell(1,0)});
    dom = struct('in',[0,1;0,1],'out',zeros(0,2));
    pvar s1_dum1 s1_dum2 s2_dum1 s2_dum2
    dpvar c1 c2
    K1m = c1*s1_dum1*s2_dum1^2 + 5*s1_dum2*s2_dum2^3 - 2*c2*s1_dum2^2*s2_dum1;
    K2m = 3*s1_dum1^2*s2_dum2 - c2*s1_dum2*s2_dum1^2 + 7*c1;
    K3m = s1_dum1 + c2*s2_dum2 + s1_dum2*s2_dum1;
    cases = {{[0,2,0,3;0,3,0,1],[0,3,0,1],{K1m;K2m}}, ...
             {[0,2,0,3;0,3,0,1;0,3,0,3],[0,3,0,1],{K1m;K2m;K3m}}, ...
             {[0,3,0,3],[0,3,0,1;0,2,0,3],{K1m,K2m}}, ...
             {[0,2,0,3],[0,2,0,3;0,3,0,1],{K1m,K2m}}};
    pvar s t s1_1 s2_1 s1_2 s2_2
    vx = s^2*t + 3*t^3 - s*t^2 + 1;
    xval = subs(vx,[s;t],[s1_1;s2_1]) * subs(vx,[s;t],[s1_2;s2_2]);
    for j=1:numel(cases)
        Kop = intop(cases{j}{3},cases{j}(1:2),vars,dom,false);
        K1 = combine_terms(Kop,2,1,true(1,2));
        check_result(Kop,K1,xval,tol,sprintf("Part 1, case %d",j));
    end

    %% Part 2: random N in {1,2,3}, several rows per variable, dpvar and polynomial kernels
    dvn = {'c1';'c2';'c3'};
    for trial = 1:40
        rng(80000+trial);
        N = randi(3);
        d = randi([2,4]);
        if N>=2,    d = min(d,3);   end
        nterms_i = randi(5,1,N);
        if N==2,    nterms_i = min(nterms_i,4);     end
        if N==3,    nterms_i = min(nterms_i,3);     end
        vnames = cell(1,N);     intmat = cell(1,N);     dumvars = cell(1,N);
        for i=1:N
            vnames{i} = ['s',num2str(i)];
            intmat{i} = rand_ties_intmat(d,nterms_i(i),0.15);
            dv = cell(d,1);
            for j=1:d,  dv{j} = [vnames{i},'_dum',num2str(j)];  end
            dumvars{i} = polynomial(dv);
        end
        vars = struct('in',{vnames},'out',{cell(1,0)});
        dom = struct('in',repmat([0,1],N,1),'out',zeros(0,2));
        allpvars = vertcat(dumvars{:});
        nt = cellfun(@(a) size(a,1),intmat);
        params = cell([nt,1]);
        use_dp = mod(trial,2)==0;
        for lin=1:numel(params)
            if use_dp
                params{lin} = rand_dpvar([1,1],dvn,allpvars,1,6);
            else
                params{lin} = rand_poly([1,1],allpvars,1,6);
            end
        end
        Kop = intop(params,intmat,vars,dom,false);
        K1 = combine_terms(Kop,d,1,true(1,N));

        tvars_p = polynomial(cellfun(@(a) ['t',num2str(a)],num2cell(1:N)','UniformOutput',false));
        vx = rand_poly([1,1],tvars_p,2,6);
        if all(isequal(vx,0)),  vx = 1+tvars_p(1);  end
        xval = 1;
        for k=1:d
            copyvars = cellfun(@(a) [a,'_',num2str(k)],vnames(:),'UniformOutput',false);
            xval = xval*subs(vx,tvars_p,polynomial(copyvars));
        end
        msg = sprintf("Part 2, trial=%d, N=%d, d=%d",trial,N,d);
        check_result(Kop,K1,xval,tol,msg);

        % (4): the order of the terms does not matter
        pm = arrayfun(@randperm,nt,'UniformOutput',false);
        if N==1,    ps = {pm{1},1};     else,   ps = pm;    end
        imat_p = cellfun(@(a,p) a(p,:),intmat,pm,'UniformOutput',false);
        Kperm = intop(params(ps{:}),imat_p,vars,dom,false);
        K1p = combine_terms(Kperm,d,1,true(1,N));
        if ~isequal(K1p.intmat,K1.intmat) || ~isequal(K1p.Z.in,K1.Z.in)
            error("case_sum_blocks failed: %s: result depends on the order of the terms.",msg)
        end
        for k=1:numel(K1.params)
            if norm(K1p.params{k}-K1.params{k},'fro')>1e-10*max(1,norm(K1.params{k},'fro'))
                error("case_sum_blocks failed: %s: parameters depend on the order of the terms.",msg)
            end
        end
    end

    disp('case_sum_blocks passed.')
end

%% ---- helpers --------------------------------------------------------------

function tf = same_op(A,B)
% Same intmat, monomials, decision variables, and parameters
tf = isequal(A.intmat,B.intmat) && isequal(A.Z.in,B.Z.in) && isequal(A.Zd,B.Zd) ...
    && isequal(size(A.params),size(B.params)) && all(cellfun(@isequal,A.params(:),B.params(:)));
end

function tf = same_kernels(A,B)
% Same intmat and the same kernels (the order of the monomials Z.in is irrelevant)
tf = isequal(A.intmat,B.intmat) && isequal(size(A.params),size(B.params));
if ~tf
    return
end
ka = kernels(A);
kb = kernels(B);
for k=1:numel(ka)
    e = cleanpoly(ka{k}-kb{k},1e-10);
    if ~isdouble(e) || any(any(abs(double(e))>1e-9))
        tf = false;
        return
    end
end
end

function tf = is_cell_row(row)
% The variables that are not multipliers form a chain, each integrated from a up
% to the next, the largest over [a,b], and each multiplier is evaluated at one of them
d = numel(row)/2;
L = row(1:2:end);       U = row(2:2:end);
is_dlt = L==U;
fre = find(~is_dlt);
up = U(~is_dlt);
tf = all(L(~is_dlt)==0) && all(ismember(L(is_dlt),fre)) && nnz(up==d+1)==1 ...
    && numel(unique(up))==numel(up) && all(ismember(up(up<=d),fre));
end

function tf = same_op_idempotent(A,B)
% A and B have the same intmat, monomials, decision variables, and parameters
tf = isequal(A.intmat,B.intmat) && isequal(A.Z.in,B.Z.in) && isequal(A.Zd,B.Zd) ...
    && isequal(size(A.params),size(B.params)) && all(cellfun(@isequal,A.params(:),B.params(:)));
end

function tf = is_canon(Kop,m)
% In each term with nonzero parameters (a combination of one intmat row per
% spatial variable), the m factors are in decreasing integration order, ties
% broken by the order along the next spatial variable
io = Kop.int_order;
N = numel(io);
nr = cellfun(@(a) size(a,1),io);
is_nz = cellfun(@(a) nnz(a)>0,Kop.params);
tf = true;
for c=1:prod(nr)
    if ~is_nz(c)
        continue
    end
    sub = cell(1,N);
    [sub{:}] = ind2sub([nr,1],c);
    key = zeros(m,N);
    for t=1:N
        % integration order of the variable each copy equals (itself, or the one a multiplier is evaluated at)
        row = Kop.intmat{t}(sub{t},:);
        is_dlt = row(1:2:end)==row(2:2:end);
        tg = 1:m;
        tg(is_dlt) = row(2*find(is_dlt)-1);
        key(:,t) = io{t}(sub{t},tg)';
    end
    for p=1:m-1
        df = key(p,:)-key(p+1,:);
        nz = find(df~=0,1);
        if ~isempty(nz) && df(nz)<0
            tf = false;
            return
        end
    end
end
end

function check_result(Kop,K1,xval,tol,msg)
% (1) same value on xval, (2) same decision variables, (3) unique sorted
% integrals and idempotent
Y0 = apply_op(Kop,xval);
Y1 = apply_op(K1,xval);
err = Y0-Y1;
if isa(err,'dpvar')
    if ~isempty(err.C) && max(abs(err.C(:)))>tol*max(1,max(abs(Y0.C(:))))
        error("case_sum_blocks failed: %s: the functional changed.",msg)
    end
else
    err = cleanpoly(err,1e-10);
    nrmref = norm(double(cleanpoly(Y0,1e-10)),'fro');
    if ~isdouble(err) || norm(double(err),'fro')>tol*max(1,nrmref)
        error("case_sum_blocks failed: %s: the functional changed.",msg)
    end
end
if ~isequal(K1.Zd,Kop.Zd)
    error("case_sum_blocks failed: %s: decision variables changed.",msg)
end
for i=1:numel(K1.intmat)
    [~,imat_u] = uniquerows_integerTable(K1.intmat{i});
    if size(imat_u,1)~=size(K1.intmat{i},1)
        error("case_sum_blocks failed: %s: integrals are not unique.",msg)
    end
end
N = numel(K1.intmat);
K2 = combine_terms(K1,max(K1.depmat.in),1,true(1,N));
if ~isequal(K2.intmat,K1.intmat) || ~isequal(K2.Z.in,K1.Z.in) || ~isequal(size(K2.params),size(K1.params))
    error("case_sum_blocks failed: %s: not unchanged by a second call.",msg)
end
for k=1:numel(K1.params)
    if ~isequal(K2.params{k},K1.params{k})
        error("case_sum_blocks failed: %s: parameters changed by a second call.",msg)
    end
end
end
