% Test_intop_uminus.m -- Tests @intop/uminus.m directly: (-Aop)(v) =
% -Aop(v), Aop+(-Aop) = 0, and -(-Aop) = Aop on a random polynomial v, all
% evaluated with apply_op.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

% Set parameters defining functional
m = 2;  n = 3;
N = 2;
dom = randi([-2,2],[N,1]) +[zeros(N,1),2*rand([N,1])];
degs = randi(3,[1,N]);
pdegs = 3;

% Set parameters defining the function
p = 2;
pdegv = 2;

% Generate a random intop object and function
Aop = rand_intop([m,n],N,dom,degs,pdegs);
pvars = Aop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n,p],pvars_full,pdegv,50);
Av = apply_op(Aop,v);

% -Aop applied to v should equal -(Aop applied to v)
Bop = -Aop;
Bv = apply_op(Bop,v);
err1 = cleanpoly(Bv - (-Av),1e-10);
nrmref = norm(full(Av.C),'fro');
if ~isdouble(err1)
    error("uminus: output is not constant.")
elseif norm(double(err1),'fro')>tol*max(1,nrmref)
    error("uminus: residual exceeds tolerance.")
end

% Aop + (-Aop) applied to v should equal 0
Cop = Aop + (-Aop);
Cv = apply_op(Cop,v);
err2 = cleanpoly(Cv,1e-10);
if ~isdouble(err2)
    error("Aop+(-Aop): output is not constant.")
elseif norm(double(err2),'fro')>tol*max(1,nrmref)
    error("Aop+(-Aop): residual exceeds tolerance (should cancel to 0).")
end

% Double negation should return the original operator's action
Dop = -(-Aop);
Dv = apply_op(Dop,v);
err3 = cleanpoly(Dv - Av,1e-10);
if ~isdouble(err3)
    error("-(-Aop): output is not constant.")
elseif norm(double(err3),'fro')>tol*max(1,nrmref)
    error("-(-Aop): residual exceeds tolerance.")
end

disp('Test_intop_uminus passed: all cases matched within tolerance.')
