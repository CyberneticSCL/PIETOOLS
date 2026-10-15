% Test_mtimes_tensopvar_N3multi.m -- Tests the intop*tensopvar composition
% when TWO OR MORE spatial variables are shared ("pass-through", N3>=2)
% between Bop's output and Bop's input (= Aop's input).
%
% Oracle: apply_op(Bop) followed by apply_op(Aop), independent of the
% composition.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Scope is Part 1 (d==1, N3 in {2,3,4}) and Part 2 (d=2, N3=2), with only
%   the shared variables. Operators with N3>=3 shared variables together
%   with extra non-shared variables are not covered.
% - Each factor's copy of a shared input variable gets a distinct dummy
%   name (_dum<j>) before the copies are combined into the input of Cop.

tol = 1e-8;

% === PART 1: d==1, ONLY the shared variables, N3 in {2,3,4} ===
nchecked = 0;
for trial = 1:20
    rng(120000+trial);

    N3 = randi([2,4]);
    shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
    vars = struct('in',{shared_vars},'out',{shared_vars});
    dom = [0,1];
    dims = randi(2,[1,2]);
    depmat = struct('in',true(1,N3),'out',true(1,N3));
    dnsty = 0.1 + 0.5*rand();

    Bop = rand_tensopvar(dims,vars,dom,depmat,1,[2,2],1,dnsty);
    pdegs = ones(1,N3);
    Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

    pvars = polynomial(Bop.vars.in(:));
    v = {rand_poly([size(Bop.ops{1,1},2),1],pvars,2,10)};
    Bv = apply_op(Bop,v);
    ABv = apply_op(Aop,Bv);
    nrmref = full(max(max(abs(double(ABv)))));
    if nrmref<=tol
        continue
    end

    Cop = mtimes_tensopvar(Aop,Bop);
    Cv = apply_op(Cop,v{1});

    err = cleanpoly(ABv-Cv,1e-10);
    if max(max(abs(err.C)))>tol*max(1,nrmref)
        error("Test_mtimes_tensopvar_N3multi failed: d==1,N3=%d case does not match ground truth (trial=%d).",N3,trial)
    end
    nchecked = nchecked+1;
end
if nchecked==0
    error("Test_mtimes_tensopvar_N3multi: Part 1 (d==1) -- every trial was skipped, no trial exercised N3>=2.")
end
fprintf('Part 1 (d==1, N3 in {2,3,4}, no extra variables) passed: %d/20 trials checked against ground truth.\n',nchecked);

% ===== PART 2: d>=2, N3=2 =======
nchecked2 = 0;
for trial = 1:20
    rng(130000+trial);

    N3 = 2; d = 2;
    vars = struct('in',{{'sShared1','sShared2'}},'out',{{'sShared1','sShared2'}});
    dom = [0,1];
    dims = randi(2,[d,2]);
    depmat = struct('in',true(d,N3),'out',true(d,N3));

    Bop = rand_tensopvar(dims,vars,dom,depmat,1,[2,2],1,0.9);
    pdegs = ones(1,N3);
    Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

    pvars = polynomial(Bop.vars.in(:));
    v = cell(1,d);
    for j=1:d
        v{j} = rand_poly([size(Bop.ops{1,j},2),1],pvars,2,10);
    end
    Bv = apply_op(Bop,v);
    ABv = apply_op(Aop,Bv);
    nrmref = full(max(max(abs(double(ABv)))));
    if nrmref<=tol
        continue
    end

    Cop = mtimes_tensopvar(Aop,Bop);
    v_dum = cell(1,d); degs_j = zeros(1,N3); v_full = 1;
    for j=1:d
        degs_j = degs_j + 1;
        pvars_new = pvars;
        for k=1:N3
            varname_k = pvar2varname(pvars(k));
            pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
        end
        v_dum{j} = subs(v{j},pvars,pvars_new);
        v_full = kron(v_full,v_dum{j});
    end
    Cv = apply_op(Cop,v_full);
    err = cleanpoly(ABv-Cv,1e-10);
    if max(max(abs(err.C)))>tol*max(1,nrmref)
        error("Test_mtimes_tensopvar_N3multi failed: d>=2,N3=2 case does not match ground truth (trial=%d).",trial)
    end
    nchecked2 = nchecked2+1;
end
if nchecked2==0
    error("Test_mtimes_tensopvar_N3multi: Part 2 (d>=2,N3=2) -- every trial was skipped, no trial exercised a nonzero case.")
end
fprintf('Part 2 (d>=2, N3=2) passed: %d/20 trials checked against ground truth.\n',nchecked2);
