% Test_intop_otimes.m -- Tests otimes for 'intop', for type [1,2]
% (elementwise along rows) and [2,2] (Kronecker product).
%
% For random functionals Aop, Bop (sharing a random subset of their
% variables), otimes(Aop,Bop,type) applied to v3 is compared to Av.*Bv
% (type [1,2]) or kron(Av,Bv) (type [2,2]), with Av = Aop(v1) and
% Bv = Bop(v2), all evaluated with apply_op. v3 = kron(v1,v2) with the
% dummy variables of the shared variables in v2 shifted past those of v1.
%
% DJ (CC), 10/15/2026: Initial coding;

% Set parameters defining functional
m1 = 2;  n1 = 2;      % row and column dimensions of functional     
m2 = 2;  n2 = 3;
N1 = 3;              % number of variables
N2 = 2;              % number of variables
dom1 = randi([-2,2],[N1,1]) +[zeros(N1,1),rand([N1,1])];
dom2 = randi([-2,2],[N2,1]) +[zeros(N2,1),rand([N2,1])];
degs1 = randi(3,[1,N1]);
degs2 = randi(3,[1,N2]);
pdegs = 3;
type = [randi(2),2];    % type of tensor product to compute (type(2)==2 required)
if type(1)==1
    m2 = m1;
end
if type(2)==1
    n2 = n1;
end

% Declare the variable names
vars1 = cell(1,N1);
for i=1:N1
    vars1{i} = ['x',num2str(i)];
end
vars2 = cell(1,N2);
n_shared = randi([0,N2]);
share_idcs2 = unique(randi(N2,[1,n_shared]));
[share_idcs1,rtn_idx] = unique(randi(N1,[1,numel(share_idcs2)]));
share_idcs2 = share_idcs2(rtn_idx);
vars2(share_idcs2) = vars1(share_idcs1);
dom2(share_idcs2,:) = dom1(share_idcs1,:);
for i=setdiff(1:N2,share_idcs2)
    vars2{i} = ['y',num2str(i)];
end

% Set parameters defining the function
p1 = 1;              % column dimension of function
p2 = 1;
pdegv = 2;

% Generate random intop objects
Aop = rand_intop([m1,n1],vars1,dom1,degs1,pdegs);
Bop = rand_intop([m2,n2],vars2,dom2,degs2,pdegs);

% Compute the product
Cop = otimes(Aop,Bop,type);

% Generate random functions on which the operators act
pvars1 = Aop.pvars.in;
pvars_full1 = polynomial(zeros(0,1));
for i=1:N1
    pvars_full1 = [pvars_full1; pvars1{i}(:)];
end
v1 = rand_poly([n1,p1],pvars_full1,pdegv,50);
pvars2 = Bop.pvars.in;
pvars_full2 = polynomial(zeros(0,1));
for i=1:N2
    pvars_full2 = [pvars_full2; pvars2{i}(:)];
end
v2 = rand_poly([n2,p2],pvars_full2,pdegv,50);

% Compute the tensor product of the functions
% Old and new variables as cellstr: subs needs a No-by-1 cellstr or polynomial array
pvarsB_old = cell(0,1);
pvarsB_new = cell(0,1);
for i=1:numel(share_idcs2)
    % Increase degree of variables in v2
    idx = share_idcs2(i);
    dA = degs1(share_idcs1(i));
    dB = degs2(idx);
    for j=1:dB
        pvarsB_old = [pvarsB_old; {[vars2{idx},'_dum',num2str(j)]}];
        pvarsB_new = [pvarsB_new; {[vars2{idx},'_dum',num2str(dA+j)]}];
    end
end
% 'new' must be an array of polynomials/doubles, not cellstr
v3 = kron(v1,subs(v2,pvarsB_old,polynomial(pvarsB_new)));

% Apply the functionals
Av = apply_op(Aop,v1);
Bv = apply_op(Bop,v2);

% Apply the new functional
Cv = apply_op(Cop,v3);

% Compare: residual must be constant and, relative to |ref_C|, below tol
tol = 1e-8;
if all(type==[1,2])
    ref_C = Av.*Bv;
elseif all(type==[2,2])
    ref_C = kron(Av,Bv);
else
    error("Test-case not supported.")
end
err_C = cleanpoly(Cv - ref_C,1e-10);
nrmref = norm(full(ref_C.C),'fro');
if ~isdouble(err_C)
    error("At least one of the outputs is not constant.")
elseif norm(double(err_C),'fro')>tol*max(1,nrmref)
    error("Outputs do not match: residual exceeds tolerance.")
else
    disp('err_C = ')
    disp(double(err_C))
end
