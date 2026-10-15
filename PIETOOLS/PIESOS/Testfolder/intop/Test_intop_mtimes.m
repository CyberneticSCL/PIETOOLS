% Test_intop_mtimes.m -- Tests left multiplication of an 'intop' by a
% constant matrix.
%
% For a random constant matrix Amat and functional Bop, (Amat*Bop)(v) is
% compared to Amat*(Bop(v)) on a random polynomial v, both evaluated with
% apply_op.
%
% DJ (CC), 10/15/2026: Initial coding;

% Set parameters defining functional
m1 = 2;  n1 = 2;      % row and column dimensions of functional     
m2 = n1;  n2 = 3;
N1 = 2;              % number of variables
N2 = 2;              % number of variables
dom = randi([-2,2],[N1,1]) +[zeros(N1,1),rand([N1,1])];
dom2 = randi([-2,2],[N2,1]) +[zeros(N2,1),rand([N2,1])];
degs = randi(3,[1,N1]);
pdegs = 3;


% Declare the variable names
vars = cell(1,N1);
for i=1:N1
    vars{i} = ['x',num2str(i)];
end

% Set parameters defining the function
p1 = 1;              % column dimension of function
p2 = 1;
pdegv = 2;

% Generate random intop objects
Amat = rand(m1,n1);
Bop = rand_intop([m2,n2],vars,dom,degs,pdegs);

% Compute the product
Cop = Amat*Bop;

% Generate random functions on which the operators act
pvars = Bop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N2
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n2,p2],pvars_full,pdegv,50);

% Apply the functionals
Bv = apply_op(Bop,v);
ABv = Amat*Bv;

% Apply the new functional
Cv = apply_op(Cop,v);

% Compare: residual must be constant and, relative to |ABv|, below tol
tol = 1e-8;
err_C = cleanpoly(Cv - ABv,1e-10);
nrmref = norm(full(ABv.C),'fro');
if ~isdouble(err_C)
    error("At least one of the outputs is not constant.")
elseif norm(double(err_C),'fro')>tol*max(1,nrmref)
    error("Outputs do not match: residual exceeds tolerance.")
else
    disp('err_C = ')
    disp(double(err_C))
end
