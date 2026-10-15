% Test_mtimes_tensopvar_zero_term.m -- Tests that composing an 'intop'
% with a 'tensopvar' in which ONE term (of several) has a zero factor
% returns the correct, generally NONZERO result (the other terms still
% contribute), rather than treating the entire operator as zero.
%
% Also covers Bop ENTIRELY zero (early return in mtimes_tensopvar,
% bypassing get_params): the result must be the zero functional of the
% correct size. Oracle: apply_op(Bop) followed by apply_op(Aop),
% independent of mtimes_tensopvar/get_params.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

for trial = 1:8
    rng(9900+trial);

    dims = [1,1; 1,1];
    vars = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
    dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
    depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
    Bop = rand_tensopvar(dims,vars,dom,depmat,3,[2,2],1,0.5);
    Bop.ops{2,1} = 0*Bop.ops{2,1};   % zero out one term's first factor

    nr = 1;
    pdegs = ones(1,numel(Bop.vars.out));
    Aop = rand_intop([nr,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

    pvars = polynomial(Bop.vars.in(:));
    v = {rand_poly([1,1],pvars(1),3,10), rand_poly([1,1],pvars(2),3,10)};
    Bv = apply_op(Bop,v);
    ABv = apply_op(Aop,Bv);
    nrmref = full(max(max(abs(double(ABv)))));
    if nrmref<=tol
        % Degenerate draw (remaining terms cancel): not informative
        continue
    end

    Cop = mtimes_tensopvar(Aop,Bop);
    v_full = kron(v{1},v{2});
    Cv = apply_op(Cop,v_full);

    err = cleanpoly(ABv-Cv,1e-10);
    if max(max(abs(err.C)))>tol*max(1,nrmref)
        error("Test_mtimes_tensopvar_zero_term failed: result does not match ground truth (trial=%d).",trial)
    end
end

% --- Complementary case: Bop entirely zero -- top-level short-circuit ---
rng(9800);
dims = [1,1; 1,1];
vars = struct('in',{{'s1','s2'}},'out',{{'s3','s4'}});
dom = struct('in',repmat([0,1],2,1),'out',repmat([0,1],2,1));
depmat = struct('in',logical(eye(2)),'out',logical(eye(2)));
Bop0 = rand_tensopvar(dims,vars,dom,depmat,3,[2,2],1,0.5);
for i=1:3
    Bop0.ops{i,1} = 0*Bop0.ops{i,1};
end
pdegs0 = ones(1,numel(Bop0.vars.out));
Aop0 = rand_intop([1,Bop0.matdim(1)],Bop0.vars.out,Bop0.dom.out,Bop0.deps.out,pdegs0);
Cop0 = mtimes_tensopvar(Aop0,Bop0);
if ~isequal(size(Cop0),[size(Aop0,1),size(Bop0,2)])
    error("Test_mtimes_tensopvar_zero_term failed: wrong size for fully-zero Bop (got %s).",mat2str(size(Cop0)))
end
pvars0 = polynomial(Bop0.vars.in(:));
v0 = {rand_poly([1,1],pvars0(1),3,10), rand_poly([1,1],pvars0(2),3,10)};
Cv0 = apply_op(Cop0,kron(v0{1},v0{2}));
if ~isdouble(Cv0) || any(any(double(Cv0)~=0))
    error("Test_mtimes_tensopvar_zero_term failed: fully-zero Bop did not yield the zero functional.")
end
fprintf('Fully-zero Bop case: passed.\n');

disp('Test_mtimes_tensopvar_zero_term passed.')
