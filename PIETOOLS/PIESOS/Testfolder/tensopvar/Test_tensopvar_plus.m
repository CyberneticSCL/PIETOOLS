% Test_tensopvar_plus.m -- tests @tensopvar/plus.m. The oracle is
% apply_op(A+B,v) against apply_op(A,v)+apply_op(B,v) on random inputs v.
% otype(1) alternates over {1,2} (pointwise vs Kronecker rows), which
% take different paths in plus.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - otype(2)=2 only: apply_op requires it, so otype(2)==1 cannot be
%   verified through this route.

ntrials = 20;
for trial=1:ntrials
% Set parameters defining tensor-PI operators
nterms = randi(3);      % number of terms in tensor-PI operator
d = randi(5);           % degree of the tensor-PI operator
N = randi(3);           % number of input variables
M = randi(3);           % number of output variables
if mod(trial,2)==1
    type = [1,2];       % pointwise rows, Kronecker columns
else
    type = [2,2];       % Kronecker rows and columns
end
dnsty = 0.25;
% Set parameters defining the function
p = 2;              % column dimension of function
pdegv = 2;

% Matrix dimensions of each factor in the tensor-PI operators
dims = randi(4,[d,2]);
if type(1)==1
    m = dims(1,1);
    dims(:,1) = m;
else
    m = prod(dims(:,1));
end
if type(2)==1
    n = dims(1,2);
    dims(:,2) = n;
else
    n = prod(dims(:,2));
end
% Input and output variables
vars = struct();
% disjoint 'x'/'y' name pools: input and output variables must not share a name (rand_sopvar rejects one name with two domains)
vars.in = arrayfun(@(i)sprintf('x%d',i),1:N,'UniformOutput',false);
vars.out = arrayfun(@(i)sprintf('y%d',i),1:M,'UniformOutput',false);
% Domains
dom = struct();
dom.in = randi([-2,2],[N,1]) +[zeros(N,1),rand([N,1])];
dom.out = randi([-2,2],[M,1]) +[zeros(M,1),rand([M,1])];
% Monomial degrees
degs = struct();
degs.in = randi(3,[1,N]);     degs.out = randi(3,[1,M]);
% Dependency arrays
depmat = struct();
depmat.in = logical(randi([0,1],d,N));
depmat.out = logical(randi([0,1],d,M));




% Generate the random tensor-PI operators and function
Aop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);
Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,type,degs,dnsty);

v = cell(1,d);
for j=1:d
    pvars_j = polynomial(Aop.ops{1,j}.vars.in(:));
    n_j = size(Aop.ops{1,j},2);
    v{j} = rand_poly([n_j,p],pvars_j,pdegv,50);
end

Cop = Aop + Bop;

Av = apply_op(Aop,v);
Bv = apply_op(Bop,v);
Cv = apply_op(Cop,v);

tol = 1e-8;
err_C = cleanpoly(Cv - (Av+Bv),1e-10);
disp('err_C = ')
disp(max(max(abs(err_C.C))))
if max(max(abs(err_C.C)))>tol
    error("Trial %d (type=[%d,%d]): outputs do not match: residual exceeds tolerance.",trial,type(1),type(2))
end
end

fprintf('Test_tensopvar_plus passed (%d trials, otype(1) in {1,2} alternating).\n',ntrials)
