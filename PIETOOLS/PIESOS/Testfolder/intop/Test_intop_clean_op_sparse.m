% Test_intop_clean_op_sparse.m -- Tests @intop/clean_op.m on LARGE SPARSE
% kernels: clean_op must work on the nonzeros of large sparse parameters,
% zero the entries below ztol (also when no parameter is entirely zero) and
% keep the others. Oracle: explicit expected sparse matrices, not clean_op
% itself.
%
% DJ (CC), 10/15/2026: Initial coding;

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
    error("Test_intop_clean_op_sparse failed: identity kernel flagged as zero.")
elseif ~isequal(Cop.params{1},speye(m))
    error("Test_intop_clean_op_sparse failed: identity kernel was changed by clean_op.")
elseif t1>5
    error("Test_intop_clean_op_sparse failed: clean_op took %.1f s for a sparse kernel with 2e4 rows (nonzero-only cleaning expected, O(nnz)).",t1)
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
    error("Test_intop_clean_op_sparse failed: sparse kernel with nonzero entries above ztol flagged as zero.")
elseif ~isequal(Cop.params{1},Pexp)
    error("Test_intop_clean_op_sparse failed: sub-tolerance entries of a sparse kernel not zeroed, or entries above ztol changed.")
end
% Custom tolerance
Cop = clean_op(Kop,0.6);
if ~isequal(Cop.params{1},sparse(ii(vv>0.6),jj(vv>0.6),vv(vv>0.6),m,m))
    error("Test_intop_clean_op_sparse failed: custom tolerance not applied to a sparse kernel.")
end

%% Case 3: every entry below the tolerance -- zero functional, dimensions kept
Kop = intop({1e-14*speye(m)},{[0,2]},{'s'},[0,1],true);
[Cop,is_zero_op] = clean_op(Kop);
if ~is_zero_op
    error("Test_intop_clean_op_sparse failed: all-below-tolerance sparse kernel not flagged as zero.")
elseif nnz(Cop.params{1})~=0 || ~isequal(size(Cop.params{1}),[m,m])
    error("Test_intop_clean_op_sparse failed: all-below-tolerance sparse kernel must be zero and keep its size.")
end

%% Case 4: dpvar kernel (coefficient matrix with a row per decision variable)
rng(150001);
q = 300;
Aop = rand_intop([2,3],{'s1'},[0,1],1,2,q);
[~,is_zero_op] = clean_op(Aop);
if is_zero_op
    error("Test_intop_clean_op_sparse failed: random dpvar functional flagged as zero.")
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
    error("Test_intop_clean_op_sparse failed: dpvar functional with some sub-tolerance coefficients flagged as zero or terms removed.")
end
for k=1:numel(Cop.params)
    if ~isequal(Cop.params{k},Cexp{k})
        error("Test_intop_clean_op_sparse failed: dpvar kernel %d: sub-tolerance coefficients not zeroed, or others changed.",k)
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
    error("Test_intop_clean_op_sparse failed: dpvar functional with all coefficients below ztol not flagged as zero.")
end
for k=1:numel(Cop.params)
    if nnz(Cop.params{k})~=0 || ~isequal(size(Cop.params{k}),size(Aop.params{k}))
        error("Test_intop_clean_op_sparse failed: dpvar kernel %d must be zero and keep its size.",k)
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
    error("Test_intop_clean_op_sparse failed: polynomial functional with some sub-tolerance coefficients flagged as zero or terms removed.")
end
for k=1:numel(Cop.params)
    if ~isequal(Cop.params{k},Cexp{k})
        error("Test_intop_clean_op_sparse failed: polynomial kernel %d: sub-tolerance coefficients not zeroed, or others changed.",k)
    end
end

disp('Test_intop_clean_op_sparse passed: large sparse kernels cleaned through their nonzeros.')
