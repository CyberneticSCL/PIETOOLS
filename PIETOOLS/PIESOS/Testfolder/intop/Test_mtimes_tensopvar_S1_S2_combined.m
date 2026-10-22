% Test_mtimes_tensopvar_S1_S2_combined.m -- Tests the intop*tensopvar
% composition when Bop has BOTH an "S2-type" output variable (new, no       % DJ (CC), 10/21/2026
% matching entry in Bop.vars.in) AND an "S1-type" input variable (consumed  % DJ (CC), 10/21/2026
% internally, no matching entry in Bop.vars.out) at once, together with
% N3>=2 shared/pass-through variables.
%
% Oracle: apply_op(Bop) followed by apply_op(Aop), independent of the
% composition.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/21/2026: Names of the variable groups S1 (input only) and S2
%   (output only) swapped, to match 'sopvar' and mtimes_tensopvar;
%
% NOTES
% - Cases with N3>=3 and d>=2 are not covered, with or without extra
%   variables.

tol = 1e-8;

% --- Part 1: d==1, N1=1,N2=1, N3 in {2,3,4,5} -- broad ---
nchecked1 = 0;
ntrials1 = 0;
for N3 = [2,3,4,5]
    for seed_ofst = [0,100,200]
        ntrials1 = ntrials1+1;
        rng(60000+N3+seed_ofst);

        d = 1;
        shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
        vars = struct('in',{[shared_vars,{'sExtraIn'}]},'out',{[shared_vars,{'sExtraOut'}]});
        depmat = struct('in',true(d,N3+1),'out',true(d,N3+1));

        Bop = rand_tensopvar([1,1],vars,[0,1],depmat,1,[2,2],1,0.5);
        pdegs = ones(1,numel(Bop.vars.out));
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
        errval = full(max(max(abs(err.C))));
        if isempty(errval); errval = 0; end
        if errval>tol*max(1,nrmref)
            error("Test_mtimes_tensopvar_S1_S2_combined failed: d==1,N3=%d case does not match ground truth (seed_ofst=%d).",N3,seed_ofst)
        end
        nchecked1 = nchecked1+1;
    end
end
if nchecked1==0
    error("Test_mtimes_tensopvar_S1_S2_combined: Part 1 -- every trial was skipped, nothing exercised.")
end
fprintf('Part 1 (d==1, N1=1,N2=1, N3 in {2,3,4,5}) passed: %d/%d trials checked against ground truth.\n',nchecked1,ntrials1);

% --- Part 2: d==2, N1=1,N2=1, N3=2 ONLY -- do not raise N3 here, see NOTES ---
nchecked2 = 0;
ntrials2 = 0;
for seed_ofst = [0,1,2]
    ntrials2 = ntrials2+1;
    N3 = 2; d = 2;
    rng(70000+N3*10+seed_ofst);

    shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
    vars = struct('in',{[shared_vars,{'sExtraIn'}]},'out',{[shared_vars,{'sExtraOut'}]});
    depmat = struct('in',true(d,N3+1),'out',true(d,N3+1));

    Bop = rand_tensopvar(repmat([1,1],d,1),vars,[0,1],depmat,1,[2,2],1,0.5);
    pdegs = ones(1,numel(Bop.vars.out));
    Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

    pvars = polynomial(Bop.vars.in(:));
    Nin_tot = numel(Bop.vars.in);
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

    % Per-factor dummy names for variables shared across >=2 factors (required for d>=2, dense depmat)
    degs_j = zeros(1,Nin_tot); v_full = 1;
    for j=1:d
        degs_j = degs_j+1;
        pvars_new = pvars;
        for k=1:Nin_tot
            varname_k = pvar2varname(pvars(k));
            pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
        end
        v_dum = subs(v{j},pvars,pvars_new);
        v_full = kron(v_full,v_dum);
    end
    Cv = apply_op(Cop,v_full);

    err = cleanpoly(ABv-Cv,1e-10);
    errval = full(max(max(abs(err.C))));
    if isempty(errval); errval = 0; end
    if errval>tol*max(1,nrmref)
        error("Test_mtimes_tensopvar_S1_S2_combined failed: d==2,N3=2 case does not match ground truth (seed_ofst=%d).",seed_ofst)
    end
    nchecked2 = nchecked2+1;
end
if nchecked2==0
    error("Test_mtimes_tensopvar_S1_S2_combined: Part 2 -- every trial was skipped, nothing exercised.")
end
fprintf('Part 2 (d==2, N1=1,N2=1, N3=2) passed: %d/%d trials checked against ground truth.\n',nchecked2,ntrials2);
