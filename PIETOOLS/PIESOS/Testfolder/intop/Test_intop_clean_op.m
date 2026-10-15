% Test_intop_clean_op.m -- Tests @intop/clean_op.m. Covers: (1) when EVERY
% kernel associated with one term along a direction is exactly zero, that
% term (row of intmat{i}) is dropped without changing the action of the
% functional; (2) sub-tolerance coefficients are zeroed without removing
% terms; (3) an entirely zero functional keeps its matdim and no direction
% collapses to 0 rows. Oracle: apply_op of the functional before vs. after
% cleaning, and the expected term counts/matdim.
%
% DJ (CC), 10/15/2026: Initial coding;

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
