function Test_intop_clean_op()
% Test_intop_clean_op.m -- one local function for each case:
%   case_basic              @intop/clean_op.m. Covers: (1) when EVERY
%                           kernel associated with one term along a
%                           direction is exactly zero, that term (row of
%                           intmat{i}) is dropped without changing the
%                           action of the
%   case_sparse             @intop/clean_op.m on LARGE SPARSE kernels:
%                           clean_op must work on the nonzeros of large
%                           sparse parameters, zero the entries below ztol
%                           (also when no parameter is entirely zero) and
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Merge Test_intop_clean_op_sparse into this test, as
%   cases;

cases = {@case_basic, @case_sparse};
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
    error('Test_intop_clean_op: %d of %d cases failed.',nfail,numel(cases));
end
fprintf('Test_intop_clean_op: all %d cases passed.\n',numel(cases));
end

%% --------------------------------------------------------------------------
function case_basic()
% case_basic -- Tests @intop/clean_op.m. Covers: (1) when EVERY
% kernel associated with one term along a direction is exactly zero, that
% term (row of intmat{i}) is dropped without changing the action of the
% functional; (2) sub-tolerance coefficients are zeroed without removing
% terms; (3) an entirely zero functional keeps its matdim and no direction
% collapses to 0 rows. Oracle: apply_op of the functional before vs. after
% cleaning, and the expected term counts/matdim.
%

    tol = 1e-8;

    % Set parameters defining functional
    m = 2;  n = 3;
    N = 2;
    dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
    pdegs = 3;
    p = 2;
    pdegv = 2;

    % Redraw degs (>=2) and the operator until some direction has >=2 terms
    % (degree 1 can collapse to 1 term in remove_badints/split_ints)
    for attempt=1:20
        degs = randi([2,3],[1,N]);
        Aop = rand_intop([m,n],N,dom,degs,pdegs);
        [ntrms_best,idir] = max(cellfun(@(a) size(a,1), Aop.intmat));
        if ntrms_best>=2
            break
        end
    end
    if ntrms_best<2
        error("Test setup error: could not find a direction with >=2 terms after 20 attempts.")
    end
    pvars = Aop.pvars.in;
    pvars_full = polynomial(zeros(0,1));
    for i=1:N
        pvars_full = [pvars_full; pvars{i}(:)];
    end
    v = rand_poly([n,p],pvars_full,pdegv,50);

    %%%% Case 1: zero out every kernel associated with term 1 along direction
    %%%% idir -- clean_op should drop that row of intmat{idir} entirely
    ntrms1 = size(Aop.intmat{idir},1);
    Bop = Aop;
    full_idcs = cellfun(@(a) 1:size(a,1), Aop.intmat, 'UniformOutput', false);
    full_idcs{idir} = 1;
    params_slice = Bop.params(full_idcs{:});
    for k=1:numel(params_slice)
        params_slice{k} = 0*params_slice{k};
    end
    Bop.params(full_idcs{:}) = params_slice;
    Bop.chkval = true;

    Cop = clean_op(Bop);
    if size(Cop.intmat{idir},1) ~= ntrms1-1
        error("clean_op: expected term 1 along direction %d to be removed (size %d -> %d), got %d.", ...
            idir, ntrms1, ntrms1-1, size(Cop.intmat{idir},1))
    end
    for i=setdiff(1:N,idir)
        if size(Cop.intmat{i},1) ~= size(Bop.intmat{i},1)
            error("clean_op: direction %d should be untouched by zeroing terms along direction %d.",i,idir)
        end
    end

    % Removing an identically-zero term must not change the action; compare with
    % the zeroed Bop, not the original Aop
    Bv = apply_op(Bop,v);
    Cv = apply_op(Cop,v);
    err1 = cleanpoly(Cv - Bv,1e-10);
    if ~isdouble(err1)
        error("clean_op: output is not constant after row removal.")
    elseif norm(double(err1),'fro')>tol*max(1,norm(double(Bv),'fro'))
        error("clean_op: functional's action changed after removing an all-zero term.")
    end

    %%%% Case 2: small (sub-tolerance) coefficients are zeroed, but no terms are
    %%%% removed if the kernel is not entirely below tolerance.
    Dop = Aop;
    ztol = 1e-6;
    noise = ztol/10;
    param1 = Dop.params{1};
    if isa(param1,'double')
        param1(1,1) = param1(1,1) + noise;
    else
        param1.C(1,1) = param1.C(1,1) + noise;
    end
    Dop.params{1} = param1;
    Eop = clean_op(Dop,ztol);
    for i=1:N
        if ~isequal(size(Eop.intmat{i}),size(Dop.intmat{i}))
            error("clean_op: term count along direction %d changed even though no kernel is entirely below tolerance.",i)
        end
    end

    %%%% Case 3: EVERY kernel is zero -- matdim must be preserved, not
    %%%% collapsed to [0,0]
    Fop = Aop;
    for k=1:numel(Fop.params)
        Fop.params{k} = 0*Fop.params{k};
    end
    Gop = clean_op(Fop);
    if ~isequal(Gop.matdim,Aop.matdim)
        error("clean_op: matdim changed (%s -> %s) when cleaning an entirely-zero functional.", ...
            mat2str(Aop.matdim),mat2str(Gop.matdim))
    end
    for i=1:N
        if isempty(Gop.intmat{i})
            error("clean_op: direction %d collapsed to 0 rows for an entirely-zero functional.",i)
        end
    end
    Gv = apply_op(Gop,v);
    if ~isdouble(Gv) || any(any(double(Gv)~=0))
        error("clean_op: entirely-zero functional did not act as the zero map after cleaning.")
    end

    disp('Test_intop_clean_op passed: zero terms removed and small coefficients cleaned without altering the functional.')
end

%% --------------------------------------------------------------------------
function case_sparse()
% case_sparse -- Tests @intop/clean_op.m on LARGE SPARSE
% kernels: clean_op must work on the nonzeros of large sparse parameters,
% zero the entries below ztol (also when no parameter is entirely zero) and
% keep the others. Oracle: explicit expected sparse matrices, not clean_op
% itself.
%

% NOTES
% - Case 1 bounds the run time (5 s at m=2e4), which excludes a dense-mask
%   implementation ('abs(P)<=ztol' on sparse P). @tensopvar/int declares
%   integrals as intop({speye(m)},..).

    m = 2e4;

    %% Case 1: identity kernel of an integral on a large space, as built by @tensopvar/int
    Kop = intop({speye(m)},{[0,2]},{'s'},[0,1],true);
    t0 = tic;
    [Cop,is_zero_op] = clean_op(Kop);
    t1 = toc(t0);
    if is_zero_op
        error("case_sparse failed: identity kernel flagged as zero.")
    elseif ~isequal(Cop.params{1},speye(m))
        error("case_sparse failed: identity kernel was changed by clean_op.")
    elseif t1>5
        error("case_sparse failed: clean_op took %.1f s for a sparse kernel with 2e4 rows (nonzero-only cleaning expected, O(nnz)).",t1)
    end

    %% Case 2: some entries below the tolerance, the kernel is not zero
    [ii,jj,vv] = find(speye(m) + sparse(1:m,[2:m,1],0.5,m,m));    % diagonal + shifted
    small = false(size(vv));        small(1:7:end) = true;          % every 7th entry
    vv(small) = 1e-14*vv(small);                                    % ... below ztol
    P = sparse(ii,jj,vv,m,m);
    Pexp = sparse(ii(~small),jj(~small),vv(~small),m,m);            % expected: the rest
    Kop = intop({P},{[0,2]},{'s'},[0,1],true);
    [Cop,is_zero_op] = clean_op(Kop);
    if is_zero_op
        error("case_sparse failed: sparse kernel with nonzero entries above ztol flagged as zero.")
    elseif ~isequal(Cop.params{1},Pexp)
        error("case_sparse failed: sub-tolerance entries of a sparse kernel not zeroed, or entries above ztol changed.")
    end
    % Custom tolerance
    Cop = clean_op(Kop,0.6);
    if ~isequal(Cop.params{1},sparse(ii(vv>0.6),jj(vv>0.6),vv(vv>0.6),m,m))
        error("case_sparse failed: custom tolerance not applied to a sparse kernel.")
    end

    %% Case 3: every entry below the tolerance -- zero functional, dimensions kept
    Kop = intop({1e-14*speye(m)},{[0,2]},{'s'},[0,1],true);
    [Cop,is_zero_op] = clean_op(Kop);
    if ~is_zero_op
        error("case_sparse failed: all-below-tolerance sparse kernel not flagged as zero.")
    elseif nnz(Cop.params{1})~=0 || ~isequal(size(Cop.params{1}),[m,m])
        error("case_sparse failed: all-below-tolerance sparse kernel must be zero and keep its size.")
    end

    %% Case 4: dpvar kernel (coefficient matrix with a row per decision variable)
    rng(150001);
    q = 300;
    Aop = rand_intop([2,3],{'s1'},[0,1],1,2,q);
    [~,is_zero_op] = clean_op(Aop);
    if is_zero_op
        error("case_sparse failed: random dpvar functional flagged as zero.")
    end
    % Some coefficients below ztol, no parameter zero: those (and only those) are zeroed
    Dop = Aop;      Cexp = cell(size(Aop.params));
    for k=1:numel(Dop.params)
        Pk = Dop.params{k};
        [ii,jj,vv] = find(Pk);
        small = false(size(vv));    small(2:3:end) = true;
        vv(small) = 1e-14;
        Pk = sparse(ii,jj,vv,size(Pk,1),size(Pk,2));
        Dop.params{k} = Pk;
        Cexp{k} = sparse(ii(~small),jj(~small),vv(~small),size(Pk,1),size(Pk,2));
    end
    [Cop,is_zero_op] = clean_op(Dop);
    if is_zero_op || ~isequal(size(Cop.params),size(Dop.params))
        error("case_sparse failed: dpvar functional with some sub-tolerance coefficients flagged as zero or terms removed.")
    end
    for k=1:numel(Cop.params)
        if ~isequal(Cop.params{k},Cexp{k})
            error("case_sparse failed: dpvar kernel %d: sub-tolerance coefficients not zeroed, or others changed.",k)
        end
    end
    Bop = Aop;
    for k=1:numel(Bop.params)
        Pk = Bop.params{k};
        Pk = 1e-14*Pk;                          % every coefficient below ztol
        Bop.params{k} = Pk;
    end
    [Cop,is_zero_op] = clean_op(Bop);
    if ~is_zero_op
        error("case_sparse failed: dpvar functional with all coefficients below ztol not flagged as zero.")
    end
    for k=1:numel(Cop.params)
        if nnz(Cop.params{k})~=0 || ~isequal(size(Cop.params{k}),size(Aop.params{k}))
            error("case_sparse failed: dpvar kernel %d must be zero and keep its size.",k)
        end
    end

    %% Case 5: plain 'polynomial' kernel, same checks
    rng(150002);
    Aop = rand_intop([2,3],{'s1'},[0,1],1,2);
    Dop = Aop;      Cexp = cell(size(Aop.params));
    for k=1:numel(Dop.params)
        Pk = Dop.params{k};
        [ii,jj,vv] = find(Pk);
        small = false(size(vv));    small(2:3:end) = true;
        vv(small) = 1e-14;
        Pk = sparse(ii,jj,vv,size(Pk,1),size(Pk,2));
        Dop.params{k} = Pk;
        Cexp{k} = sparse(ii(~small),jj(~small),vv(~small),size(Pk,1),size(Pk,2));
    end
    [Cop,is_zero_op] = clean_op(Dop);
    if is_zero_op || ~isequal(size(Cop.params),size(Dop.params))
        error("case_sparse failed: polynomial functional with some sub-tolerance coefficients flagged as zero or terms removed.")
    end
    for k=1:numel(Cop.params)
        if ~isequal(Cop.params{k},Cexp{k})
            error("case_sparse failed: polynomial kernel %d: sub-tolerance coefficients not zeroed, or others changed.",k)
        end
    end

    disp('case_sparse passed: large sparse kernels cleaned through their nonzeros.')
end
