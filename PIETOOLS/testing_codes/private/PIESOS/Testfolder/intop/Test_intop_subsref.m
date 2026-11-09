% Test_intop_subsref.m -- Tests @intop/subsref.m: row and column indexing,
% bounds checks, and dot indexing.
%
% Oracle: Kop(ridcs,:)*v is compared with (Kop*v)(ridcs,:), and
% Kop(:,cidcs)*v(cidcs,:) with Kop*v for v zeroed outside cidcs, both
% evaluated with apply_op.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - n > m on purpose, with a column index exceeding m: the column bound is
%   checked against size(obj,2), not size(obj,1).

tol = 1e-8;

% Set parameters defining functional; n > m on purpose (see NOTES)
m = 2;  n = 5;
N = 2;
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
degs = randi(3,[1,N]);
pdegs = 3;
p = 2;
pdegv = 2;

Kop = rand_intop([m,n],N,dom,degs,pdegs);
pvars = Kop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n,p],pvars_full,pdegv,50);
Kv = apply_op(Kop,v);

%%%% Row indexing: Kop(ridcs,:)*v should equal (Kop*v)(ridcs,:)
ridcs = [2,1];      % also checks reordering, not just selection
Ksub_r = Kop(ridcs,:);
if ~isequal(size(Ksub_r),[numel(ridcs),n])
    error("Row indexing: size of Kop(ridcs,:) is incorrect.")
end
Kv_r = apply_op(Ksub_r,v);
Kv_rsub = Kv(ridcs,:);
err_r = cleanpoly(Kv_r - Kv_rsub,1e-10);
nrmref_r = norm(full(Kv_rsub.C),'fro');
if ~isdouble(err_r)
    error("Row indexing: output is not constant.")
elseif norm(double(err_r),'fro')>tol*max(1,nrmref_r)
    error("Row indexing: residual exceeds tolerance.")
end

%%%% Column indexing, with a column index > m (only possible because n > m)
cidcs = [1,4];
if max(cidcs)<=m
    error("Test setup error: expected a column index exceeding m to exercise the fix.")
end
Ksub_c = Kop(:,cidcs);
if ~isequal(size(Ksub_c),[m,numel(cidcs)])
    error("Column indexing: size of Kop(:,cidcs) is incorrect.")
end
% Column-restricted operator on the corresponding rows of v equals the full
% operator on v with all other rows zeroed
v_zeroed = 0*v;
v_zeroed(cidcs,:) = v(cidcs,:);
Kv_zeroed = apply_op(Kop,v_zeroed);
Kv_c = apply_op(Ksub_c,v(cidcs,:));
err_c = cleanpoly(Kv_c - Kv_zeroed,1e-10);
nrmref_c = norm(full(Kv_zeroed.C),'fro');
if ~isdouble(err_c)
    error("Column indexing: output is not constant.")
elseif norm(double(err_c),'fro')>tol*max(1,nrmref_c)
    error("Column indexing: residual exceeds tolerance.")
end

%%%% Out-of-range indices should still error
threw_error = false;
try
    Ksub_bad = Kop(:,n+1); %#ok<NASGU>
catch
    threw_error = true;
end
if ~threw_error
    error("Column indexing: out-of-range index should have thrown an error.")
end

%%%% Dot-indexing (field access) should still work as a passthrough
if ~isequal(Kop.matdim,[m,n])
    error("Dot-indexing: Kop.matdim does not match [m,n].")
end

disp('Test_intop_subsref passed: row/column indexing and bounds checks behave correctly.')
