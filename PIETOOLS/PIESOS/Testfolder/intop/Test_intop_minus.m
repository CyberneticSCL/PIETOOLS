% Test_intop_minus.m -- Tests minus for 'intop'.
%
% For random functionals Aop, Bop on N variables, (Aop-Bop)(v) is compared
% to Aop(v)-Bop(v) on a random polynomial v, both evaluated with apply_op.
%
% DJ (CC), 10/15/2026: Initial coding;

% Set parameters defining functional
m = 2;  n = 3;      % row and column dimensions of functional              
N = 3;              % number of variables
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
degs = randi(4,[1,N]);
pdegs = 3;

% Set parameters defining the function
p = 2;              % column dimension of function
pdegv = 2;

% Generate random intop objects
Aop = rand_intop([m,n],N,dom,degs,pdegs);
Bop = rand_intop([m,n],N,dom,degs,pdegs);

% Compute the difference
Cop = Aop-Bop;

% Generate random functions on which the operators act
pvars = Aop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n,p],pvars_full,pdegv,50);

% Apply the functionals
Av = apply_op(Aop,v);
Bv = apply_op(Bop,v);

% Apply the new functional
Cv = apply_op(Cop,v);

% Compare: residual must be constant and, relative to |Av-Bv|, below tol
tol = 1e-8;
AmBv = Av-Bv;
err_C = cleanpoly(Cv - AmBv,1e-10);
nrmref = norm(full(AmBv.C),'fro');
if ~isdouble(err_C)
    error("At least one of the outputs is not constant.")
elseif norm(double(err_C),'fro')>tol*max(1,nrmref)
    error("Outputs do not match: residual exceeds tolerance.")
else
    disp('err_C = ')
    disp(double(err_C))
end
