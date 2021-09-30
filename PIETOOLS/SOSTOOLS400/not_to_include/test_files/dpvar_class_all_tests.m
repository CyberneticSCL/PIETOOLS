function logval = dpvar_class_all_tests()
% test creation, create random dpvar
test_dpvar_fidelity_DPcommon_basis;
test_dpvar_fidelity_DPcompress;
test_dpvar_fidelity_dpvar2poly;
test_dpvar_fidelity_poly2dpvar;
test_dpvar_fidelity_mtimes;
test_dpvar_fidelity_plusminus;
test_dpvar_fidelity_transpose;
test_dpvar_fidelity_horzcat;
test_dpvar_fidelity_vertcat;
test_dpvar_fidelity_blkdiag;
test_dpvar_fidelity_subsref;
test_dpvar_fidelity_subsasgn;
test_dpvar_fidelity_int;
test_dpvar_fidelity_subs;
test_dpvar_fidelity_varswap;
test_dpvar_fidelity_jacobian;
test_dpvar_fidelity_rdivide;

end

function testDp = random_dpvar(varargs)
if nargin==0
    m = randi([1,20],1); n = randi([1,20],1);
    ndegmat = randi([0,10],1);
    npvars = randi([0,10],1);
    ndvars = randi([0,5000],1);
    sparsity = rand(1);
elseif nargin==6
    m = varargs{1}; n = varargs{2}; ndegmat = varargs{3};
    npvars = varargs{4}; ndvars = varargs{5}; sparsity = varargs{6};
else
    error("Input to random_dpvar should be empty or 6 inputs corresponding to rows, cols, degmat_rows, pvar numbers, dvar numbers and sparsity value between (0,1)");
end

dvarname = {}; varname = {};
for i=1:npvars
    varname = [varname, ['pvar_',num2str(i)]];
end
for i=1:ndvars
    dvarname = [dvarname, ['coeff_',num2str(i)]];
end

[~,idxp] = sort(rand(1,npvars));
[~,idxd] = sort(rand(1,ndvars));

matdim = [m,n];
degmat = floor(5*sprand(ndegmat,npvars,sparsity));

C = sprand((ndvars+1)*m,ndegmat*n,sparsity);

testDp = dpvar(C,degmat,varname(idxp),dvarname(idxd),matdim);
end