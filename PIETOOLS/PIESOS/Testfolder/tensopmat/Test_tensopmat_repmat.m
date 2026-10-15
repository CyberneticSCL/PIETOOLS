% Test_tensopmat_repmat.m -- tests @tensopmat/repmat.m.
%
% B = repmat(A,r1,r2) tiles the WHOLE p x d block matrix r1 x r2 times, so
% B.ops = repmat(A.ops,r1,r2) (MATLAB's own cell repmat). Feeding B the
% input X tiled r2 times column-wise (X_tiled = repmat(X,1,r2)), each of
% the r1 row-tiles sees the same d column-blocks repeated r2 times, so
% every row-tile evaluates to r2*apply_op(A,X):
% Y_tiled = r2*repmat(Y,r1,1).

% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

p = randi([1,3]);  d = randi([1,3]);
N = randi([2,3]);  M = randi([2,3]);
dom = [0,1];

[Aop,info] = rand_tensopmat(p,d,N,M,dom,0.7);
X = cell(1,d);
for j=1:d
    X{j} = rand_poly([info.dims_in(j),1],info.vars_in_pool(info.sel_in{j}),1,5);
end
Y = apply_op(Aop,X);

r1 = randi([2,3]);  r2 = randi([2,3]);
Bop = repmat(Aop,r1,r2);

if ~isequal(size(Bop),[size(Aop,1)*r1,size(Aop,2)*r2])
    error("Test_tensopmat_repmat failed: size(repmat(A,r1,r2)) mismatch.")
end
if ~isequal(Bop.dims.out,repmat(Aop.dims.out,r1,1)) || ~isequal(Bop.dims.in,repmat(Aop.dims.in,r2,1))
    error("Test_tensopmat_repmat failed: dims.out/dims.in not tiled correctly.")
end
if ~isequal(Bop.depmat.out,repmat(Aop.depmat.out,r1,1)) || ~isequal(Bop.depmat.in,repmat(Aop.depmat.in,r2,1))
    error("Test_tensopmat_repmat failed: depmat.out/depmat.in not tiled correctly.")
end

Xtiled = repmat(X,1,r2);
Ytiled = apply_op(Bop,Xtiled);
err = cleanpoly(Ytiled-r2*repmat(Y,r1,1),1e-10);
if max(max(abs(err.C)))>tol
    error("Test_tensopmat_repmat failed: apply(repmat(A,r1,r2),Xtiled) does not match r2*repmat(apply(A,X),r1,1).")
end

% r1=1,r2=1 must be a no-op
Bop1 = repmat(Aop,1,1);
Y1 = apply_op(Bop1,X);
err1 = cleanpoly(Y1-Y,1e-10);
if max(max(abs(err1.C)))>tol
    error("Test_tensopmat_repmat failed: repmat(A,1,1) does not reproduce A.")
end

disp('Test_tensopmat_repmat passed.')
