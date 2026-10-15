% Test_intop_unique_intmat.m -- Tests @intop/unique_intmat.m. Every row of
% intmat{i} is duplicated along every direction (so consolidation is
% exercised along all N directions at once), and the test checks that no
% duplicate rows remain and that the value of the functional on a random
% polynomial (apply_op) is unchanged by consolidation.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

% Set parameters defining functional
m = 2;  n = 3;
N = 2;
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
degs = randi(3,[1,N]);
pdegs = 3;
p = 2;
pdegv = 2;

Aop = rand_intop([m,n],N,dom,degs,pdegs);
pvars = Aop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n,p],pvars_full,pdegv,50);

%%%% Duplicate every row of intmat{i} along every direction, and expand params to match
Aop_dup = Aop;
full_idcs = cell(1,N);
for i=1:N
    Aop_dup.intmat{i} = [Aop.intmat{i}; Aop.intmat{i}];
    full_idcs{i} = [1:size(Aop.intmat{i},1), 1:size(Aop.intmat{i},1)];
end
Aop_dup.params = Aop.params(full_idcs{:});
Aop_dup.chkval = false;

%%%% The functional's action must be unaffected by consolidating duplicate integrals
Av_dup = apply_op(Aop_dup,v);
Aop_uniq = unique_intmat(Aop_dup);
Av_uniq = apply_op(Aop_uniq,v);

for i=1:N
    if size(Aop_uniq.intmat{i},1) > size(Aop_dup.intmat{i},1)
        error("unique_intmat: number of terms along direction %d did not decrease despite duplicated rows.",i)
    elseif size(Aop_uniq.intmat{i},1) ~= size(unique(Aop_uniq.intmat{i},'rows'),1)
        error("unique_intmat: direction %d still contains duplicate rows after consolidation.",i)
    end
end

err = cleanpoly(Av_uniq - Av_dup,1e-10);
if ~isdouble(err)
    error("unique_intmat: output is not constant.")
elseif norm(double(err),'fro')>tol*max(1,norm(double(Av_dup),'fro'))
    error("unique_intmat: functional's action changed after consolidating duplicate integrals.")
end

%%%% unique_intmat on an already-unique operator leaves intmat and the action unchanged
Aop_uniq2 = unique_intmat(Aop_uniq);
for i=1:N
    if ~isequal(Aop_uniq2.intmat{i},Aop_uniq.intmat{i})
        error("unique_intmat: applying to an already-unique operator changed intmat{%d}.",i)
    end
end
Av_uniq2 = apply_op(Aop_uniq2,v);
err2 = cleanpoly(Av_uniq2 - Av_uniq,1e-10);
if ~isdouble(err2) || norm(double(err2),'fro')>tol*max(1,norm(double(Av_uniq),'fro'))
    error("unique_intmat: idempotence check failed.")
end

disp('Test_intop_unique_intmat passed: duplicate integrals consolidated without changing the functional.')
