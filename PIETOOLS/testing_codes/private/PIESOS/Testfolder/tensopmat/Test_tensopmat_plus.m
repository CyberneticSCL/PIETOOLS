% Test_tensopmat_plus.m -- tests plus/minus/uminus for 'tensopmat'
% objects, via apply_op(A+A,X) == 2*apply_op(A,X), apply_op(A-A,X) == 0,
% and apply_op(-A,X) == -apply_op(A,X).
%
% Sweeps otype(1) in {1,2}, nterms and degree per block, p/d/N/M down to
% 1x1, and a sparser 'dnsty' on some reps; 20 fresh random reps.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES
% - Checks are self-consistent (A against itself), so mixing otype(1)
%   across blocks carries none of the apply_op output renaming issue
%   described in rand_tensopmat.m / Test_tensopmat_mtimes.m; no
%   normalization is needed.

tol = 1e-8;
nreps = 20;

for rep_idx=1:nreps

p = randi([1,4]);  d = randi([1,4]);
N = randi([1,4]);  M = randi([1,4]);
dom = [0,1];
if rand()<0.3
    dnsty = 0.3+0.4*rand();
else
    dnsty = [];
end

[Cop,info] = rand_tensopmat(p,d,N,M,dom,0.7,dnsty,struct(),[],[],'rand',3,3);

X = cell(1,d);
for j=1:d
    X{j} = rand_poly([info.dims_in(j),1],info.vars_in_pool(info.sel_in{j}),1,5);
end

Y = apply_op(Cop,X);

% Test plus: A+A
Ysum = apply_op(Cop+Cop,X);
err_plus = cleanpoly(Ysum-2*Y,1e-10);
if max(max(abs(err_plus.C)))>tol
    error("Test_tensopmat_plus failed: A+A does not match 2*apply(A,X).")
end

% Test minus: A-A
Ydiff = apply_op(Cop-Cop,X);
err_minus = cleanpoly(Ydiff,1e-10);
if max(max(abs(err_minus.C)))>tol
    error("Test_tensopmat_plus failed: A-A does not evaluate to 0.")
end

% Test uminus
Yneg = apply_op(-Cop,X);
err_uminus = cleanpoly(Yneg+Y,1e-10);
if max(max(abs(err_uminus.C)))>tol
    error("Test_tensopmat_plus failed: -A does not match -apply(A,X).")
end

end

disp('Test_tensopmat_plus passed.')
