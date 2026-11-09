% Test_intop_times.m -- Tests @intop/times.m (the '.*' operator): the
% scalar-multiplier case and the elementwise same-size case, both operand
% orders (M.*Kop and Kop.*M).
%
% Oracle: Case 1 compares apply_op(a.*Bop,v) with a*apply_op(Bop,v). Case 2
% compares the kernels of C2 and Bop, evaluated from the definition of the
% class (intop_eval_kernel) at random values of the dummy variables, as
% Amat.*(kernel of Bop).
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Case 2 is not checked through apply_op: scaling the kernel entrywise by
%   a non-scalar Amat does not commute with the summation over columns.
% - Tolerances are relative: |Bv| routinely reaches 1e7-1e10 (domain
%   lengths up to ~4, pdegs=3, pdegv=2), where an absolute 1e-8 is
%   unmeetable in double precision.

tol = 1e-8;

% Set parameters defining functional
m = 2;  n = 3;      % row and column dimensions of functional
N = 2;              % number of variables
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
degs = randi(3,[1,N]);
pdegs = 3;

% Set parameters defining the function
p = 2;              % column dimension of function
pdegv = 2;

% Generate a random intop object
Bop = rand_intop([m,n],N,dom,degs,pdegs);

% Generate a random function on which the operator acts
pvars = Bop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n,p],pvars_full,pdegv,50);
Bv = apply_op(Bop,v);

%%%% Case 1: scalar multiplier, both orders
a = randn;
C1 = a.*Bop;
C1v = apply_op(C1,v);
err1 = cleanpoly(C1v - a*Bv,1e-10);
if ~isdouble(err1)
    error("Case 1 (a.*Bop): output is not constant.")
elseif norm(double(err1),'fro')>tol*max(1,norm(double(a*Bv),'fro'))
    error("Case 1 (a.*Bop): residual exceeds tolerance.")
end

C1b = Bop.*a;
C1bv = apply_op(C1b,v);
err1b = cleanpoly(C1bv - a*Bv,1e-10);
if ~isdouble(err1b)
    error("Case 1b (Bop.*a): output is not constant.")
elseif norm(double(err1b),'fro')>tol*max(1,norm(double(a*Bv),'fro'))
    error("Case 1b (Bop.*a): residual exceeds tolerance.")
end

%%%% Case 2: elementwise same-size multiplier, both orders
Amat = randn(m,n);
C2 = Amat.*Bop;
% Dummy variable values at which to evaluate the kernels
vn = Bop.varnames.in_dum;
dumnames = [cell(1,0),vn{:}];
vals = containers.Map(dumnames,num2cell(rand(1,numel(dumnames))));
for k=1:numel(Bop.params)
    Kc = intop_eval_kernel(C2,k,vals,zeros(0,1));
    Kb = intop_eval_kernel(Bop,k,vals,zeros(0,1));
    if norm(Kc-Amat.*Kb,'fro')>tol*max(1,norm(Kb,'fro'))
        error("Case 2 (Amat.*Bop): kernel %d does not match Amat.*Bop.params{%d}.",k,k)
    end
end

C2b = Bop.*Amat;
for k=1:numel(Bop.params)
    Kc = intop_eval_kernel(C2b,k,vals,zeros(0,1));
    Kb = intop_eval_kernel(Bop,k,vals,zeros(0,1));
    if norm(Kc-Amat.*Kb,'fro')>tol*max(1,norm(Kb,'fro'))
        error("Case 2b (Bop.*Amat): kernel %d does not match Amat.*Bop.params{%d}.",k,k)
    end
end

disp('Test_intop_times passed: all cases matched within tolerance.')
