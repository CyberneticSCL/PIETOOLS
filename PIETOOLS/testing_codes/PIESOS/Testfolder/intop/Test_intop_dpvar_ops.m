% Test_intop_dpvar_ops.m -- Tests 'dpvar'-kernelled 'intop' objects under
% plus, mtimes (multiplier*functional), and otimes. See also the
% dpvar-kernel check in Test_intop_mtimes_tensopvar.m.
%
% Oracle: for each operation, the decision variables are set to fixed
% random numeric values d0 either (a) before the operation, giving plain
% 'polynomial' operands that go through the q=0 code path, or (b) after the
% operation on the symbolic 'dpvar' objects. The operations are linear in
% the kernel coefficients, so (a) and (b) must agree; both are applied to
% the same random function and compared numerically.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - dpvar*dpvar is not affine and is disallowed by the class, so every
%   binary operation pairs one dpvar-kernelled intop with a plain operand.
%   'plus' is the exception: both operands are dpvar-kernelled (sharing one
%   set of q decision variables, as rand_intop guarantees).

tol = 1e-7;
nreps = 20;
q = 60;         % number of decision variables

for rep=1:nreps
rng(10000+rep);

%% ---- PLUS: two dpvar-kernelled intops (same shared q decision vars) ----
m = randi(3);   n = randi(3);   N = randi([1,5]);
dom = randi([-2,2],[N,1]) + [zeros(N,1),rand([N,1])];
degs = randi(3,[1,N]);
pdegs = randi(3);
Aop = rand_intop([m,n],N,dom,degs,pdegs,q);
Bop = rand_intop([m,n],N,dom,degs,pdegs,q);

dvn = Aop.Zd;
d0 = randn(numel(dvn),1);
d0map = containers.Map(dvn,num2cell(d0));

Aop_num = eval_intop_dpvar(Aop,d0map);
Bop_num = eval_intop_dpvar(Bop,d0map);

Cop_dpvar = Aop + Bop;
Cop_truth = Aop_num + Bop_num;

pvars = Aop.pvars.in;
pvars_full = polynomial(zeros(0,1));
for i=1:N
    pvars_full = [pvars_full; pvars{i}(:)];
end
v = rand_poly([n,1],pvars_full,2,20);

Cv_dpvar_num = eval_dpvar_at(apply_op(Cop_dpvar,v),d0map);
Cv_truth = apply_op(Cop_truth,v);
err = cleanpoly(Cv_dpvar_num - Cv_truth,1e-10);
nrmref = full(max(max(abs(double(Cv_truth)))));
if ~isdouble(err)
    error("Test_intop_dpvar_ops (plus, rep %d): output is not constant.",rep)
elseif norm(double(err),'fro')>tol*max(1,nrmref)
    error("Test_intop_dpvar_ops (plus, rep %d): residual exceeds tolerance.",rep)
end

%% ---- MTIMES: plain double multiplier * dpvar-kernelled functional ----
p = randi(3);
Amat = randn(p,m);
Cop_dpvar = Amat*Aop;
Cop_truth = Amat*Aop_num;
Cv_dpvar_num = eval_dpvar_at(apply_op(Cop_dpvar,v),d0map);
Cv_truth = apply_op(Cop_truth,v);
err = cleanpoly(Cv_dpvar_num - Cv_truth,1e-10);
nrmref = full(max(max(abs(double(Cv_truth)))));
if ~isdouble(err)
    error("Test_intop_dpvar_ops (mtimes, rep %d): output is not constant.",rep)
elseif norm(double(err),'fro')>tol*max(1,nrmref)
    error("Test_intop_dpvar_ops (mtimes, rep %d): residual exceeds tolerance.",rep)
end

%% ---- OTIMES: dpvar-kernelled Aop otimes plain-polynomial intop ----
m2 = m;          n2 = randi(2);   N2 = randi([1,4]);     % type(1)=1 needs m1==m2
vars2 = cell(1,N2);
for i=1:N2
    vars2{i} = ['y',num2str(i)];
end
dom2 = randi([-2,2],[N2,1]) + [zeros(N2,1),rand([N2,1])];
degs2 = randi(2,[1,N2]);
Bop2 = rand_intop([m2,n2],vars2,dom2,degs2,randi(2));  % q=0: plain polynomial

Cop_dpvar = otimes(Aop,Bop2,[1,2]);
Cop_truth = otimes(Aop_num,Bop2,[1,2]);

pvars2 = Bop2.pvars.in;
pvars_full2 = polynomial(zeros(0,1));
for i=1:N2
    pvars_full2 = [pvars_full2; pvars2{i}(:)];
end
v2 = rand_poly([n2,1],pvars_full2,2,20);
v3 = kron(v,v2);    % Aop and Bop2 share no variables by construction

Cv_dpvar_num = eval_dpvar_at(apply_op(Cop_dpvar,v3),d0map);
Cv_truth = apply_op(Cop_truth,v3);
err = cleanpoly(Cv_dpvar_num - Cv_truth,1e-10);
nrmref = full(max(max(abs(double(Cv_truth)))));
if ~isdouble(err)
    error("Test_intop_dpvar_ops (otimes, rep %d): output is not constant.",rep)
elseif norm(double(err),'fro')>tol*max(1,nrmref)
    error("Test_intop_dpvar_ops (otimes, rep %d): residual exceeds tolerance.",rep)
end

end

fprintf('Test_intop_dpvar_ops passed: %d reps of plus/mtimes/otimes on dpvar-kernelled intop (q=%d).\n',nreps,q)
