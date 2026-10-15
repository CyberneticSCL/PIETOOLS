% Test_common_intmat.m -- Tests common_intmat: A and B expressed in the
% same integrals must act on a random function exactly as before. Oracle:
% apply_op of the original operators.
%
% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Relative tolerance: |Av|,|Bv| reach ~1e8 (degs up to 4 in N=3
%   directions, domain length up to ~4), so an absolute tol fails on
%   roundoff.

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

% Express A and B in terms of the same integrals
[Aop2,Bop2] = common_intmat(Aop,Bop);

% Apply the new functionals
Av2 = apply_op(Aop2,v);
Bv2 = apply_op(Bop2,v);

% Compare
% Residual must be a constant that is (relatively) zero
tol = 1e-8;
err_A = cleanpoly(Av2 - Av,1e-10);
err_B = cleanpoly(Bv2 - Bv,1e-10);
if ~isdouble(err_A) || ~isdouble(err_B)
    error("The two outputs do not match.")
elseif norm(double(err_A),'fro')>tol*max(1,norm(double(Av),'fro')) ...
        || norm(double(err_B),'fro')>tol*max(1,norm(double(Bv),'fro'))
    error("The two outputs do not match: residual exceeds tolerance.")
else
    disp('err_A=')
    disp(double(err_A))
    disp('err_B=')
    disp(double(err_B))
end
