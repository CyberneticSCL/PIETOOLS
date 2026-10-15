% Test_tensopmat_repelem.m -- tests @tensopmat/repelem.m.
%
% B = repelem(A,r1,r2) is the 'tensopmat' representing kron(A,ones(r1,r2)):
% every row and column of the matrix defined by A is repeated r1, r2 times,
% in place, also within the blocks (each block keeps its variables, only
% its dimensions, and those of the object, change). So for the inputs
% Xr{j} = repelem(X{j},r2,1) (each element of X{j} repeated r2 times),
% apply_op(B,Xr) = r2*repelem(apply_op(A,X),r1,1), which is checked here,
% together with the dimensions, the dependency arrays, and the size of the
% blocks. The repeated rows and columns are taken by indexing, not by the
% builtin repelem, which errors for a 'polynomial' (SOSTOOLS400). p, d
% include 1, to cover single row-/column-blocks, and the blocks have
% several rows and columns (rand_tensopmat).

% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

for rep_idx=1:20

p = randi([1,3]);  d = randi([1,3]);
N = randi([2,3]);  M = randi([2,3]);
dom = [0,1];

[Aop,info] = rand_tensopmat(p,d,N,M,dom,0.7,[],struct(),[],[],'rand',3,3);       % pointwise and Kronecker blocks, several terms
X = cell(1,d);
for j=1:d
    X{j} = rand_poly([info.dims_in(j),1],info.vars_in_pool(info.sel_in{j}),1,5);
end
Y = apply_op(Aop,X);

r1 = randi([2,3]);  r2 = randi([2,3]);
Bop = repelem(Aop,r1,r2);

if ~isequal(size(Bop),[size(Aop,1)*r1,size(Aop,2)*r2])
    error("Test_tensopmat_repelem failed: size(repelem(A,r1,r2)) mismatch.")
end
if ~isequal(Bop.dims.out,r1*Aop.dims.out) || ~isequal(Bop.dims.in,r2*Aop.dims.in)
    error("Test_tensopmat_repelem failed: dims.out/dims.in not repeated correctly.")
end
if ~isequal(Bop.depmat.out,Aop.depmat.out) || ~isequal(Bop.depmat.in,Aop.depmat.in) || ~isvalid(Bop)
    error("Test_tensopmat_repelem failed: depmat.out/depmat.in not repeated correctly.")
end

% each input element repeated r2 times, each output element r1 times
Xr = cell(1,d);
for j=1:d
    Xr{j} = X{j}(repelem(1:size(X{j},1),r2),:);
end
Yr = apply_op(Bop,Xr);
Yexpected = r2*Y(repelem(1:size(Y,1),r1),:);
err = cleanpoly(Yr-Yexpected,1e-10);
if max(max(abs(err.C)))>tol
    error("Test_tensopmat_repelem failed: apply(repelem(A,r1,r2),Xr) does not match r2*repelem(apply(A,X),r1,1).")
end
% r1=1,r2=1 must be a no-op
Bop1 = repelem(Aop,1,1);
Y1 = apply_op(Bop1,X);
err1 = cleanpoly(Y1-Y,1e-10);
if max(max(abs(err1.C)))>tol
    error("Test_tensopmat_repelem failed: repelem(A,1,1) does not reproduce A.")
end

end

% A functional ('intop'): x -> int_0^1 x, repeated 2 x 3
pvar s
Kint = tensopmat(intop({1},{[0,2]},{'s'},[0,1],true));
Brep = repelem(Kint,2,3);
x = 1+s+s^2;        Ix = 1+1/2+1/3;
Y = polynomial(apply_op(Brep,{[x;x;x]}));
err = cleanpoly(Y-[3*Ix;3*Ix],1e-12);
if ~isequal(size(Brep),[2,3]) || ~isvalid(Brep) || max(abs(err.C(:)))>tol
    error("Test_tensopmat_repelem failed: repelem of a functional.")
end
% A block of type 'polynomial', acting on R: each row repeated twice
vn = struct('in',{cell(1,0)},'out',{{'s'}});    dn = struct('in',zeros(0,2),'out',[0,1]);
Cnat = tensopmat({[1+s,s;2,s^2]},vn,dn);
Bnat = repelem(Cnat,2,1);
Ynat = apply_op(Bnat,{[1;2]});
Yref = [1+3*s;1+3*s;2+2*s^2;2+2*s^2];
err = cleanpoly(Ynat-Yref,1e-12);
if ~isequal(size(Bnat),[4,2]) || ~isvalid(Bnat) || ~isequal(Bnat.dims.out,4) || max(abs(err.C(:)))>tol
    error("Test_tensopmat_repelem failed: repelem of a matrix does not repeat its rows.")
end
disp('Test_tensopmat_repelem passed.')
