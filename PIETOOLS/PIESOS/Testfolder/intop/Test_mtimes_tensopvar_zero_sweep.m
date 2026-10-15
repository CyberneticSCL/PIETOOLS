% Test_mtimes_tensopvar_zero_sweep.m -- Randomized sweep of the
% intop*tensopvar composition across varying numbers of spatial variables,
% factors, terms, and a wide range of coefficient densities (to generate a
% diverse mix of structural zero patterns across terms/factors), checking
% the result against an independent ground truth every time.
%
% Companion to the single hand-picked scenario of
% Test_mtimes_tensopvar_zero_term.m: no structural zero pattern may make
% mtimes_tensopvar silently return a wrong (typically all-zero) result.
% Also exercises the checks in get_params of mtimes_tensopvar.m (error ids
% 'get_params:inconsistent_zero_slot' and
% 'get_params:inconsistent_params_type'), which fail loudly if they fire in
% error. Oracle: apply_op(Bop) followed by apply_op(Aop), independent of
% mtimes_tensopvar/get_params.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Each factor's copy of a shared input variable gets its own dummy name
%   (_dum<j>) before the copies are combined (required when >=2 factors
%   depend on the same input variable).

tol = 1e-8;
ntrials = 150;
nskipped = 0;
nchecked = 0;

for trial = 1:ntrials
    rng(70000+trial);

    N = randi([1,2]);
    M = randi([1,2]);
    d = randi([1,3]);
    nterms = randi([2,4]);
    dims = randi(2,[d,2]);
    % Force exactly one variable shared between vars.in/vars.out (N3>=1).
    vars_in = [arrayfun(@(k) sprintf('sIn%d',k),1:N-1,'UniformOutput',false),{'sShared'}];
    vars_out = [arrayfun(@(k) sprintf('sOut%d',k),1:M-1,'UniformOutput',false),{'sShared'}];
    vars = struct('in',{vars_in},'out',{vars_out});
    dom = [0,1];
    depmat = struct('in',logical(randi([0,1],d,N)),'out',logical(randi([0,1],d,M)));
    depmat.in(1,~any(depmat.in,1)) = true;
    depmat.out(1,~any(depmat.out,1)) = true;
    dnsty = 0.1 + 0.5*rand();   % vary density to mix partial- and full-zero patterns

    try
        Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],1,dnsty);
    catch
        nskipped = nskipped+1;
        continue
    end
    [~,is_zeroB] = clean_op(Bop);
    if is_zeroB
        continue   % trivially-zero case already covered by Test_mtimes_tensopvar_zero_term.m
    end

    pdegs = ones(1,numel(Bop.vars.out));
    try
        Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);
    catch
        nskipped = nskipped+1;
        continue
    end

    pvars = polynomial(Bop.vars.in(:));
    v = cell(1,d);
    for j=1:d
        n_j = size(Bop.ops{1,j},2);
        is_dep_j = depmat.in(j,:);
        v{j} = rand_poly([n_j,1],pvars(is_dep_j),2,10);
    end
    Bv = apply_op(Bop,v);
    ABv = apply_op(Aop,Bv);
    nrmref = full(max(max(abs(double(ABv)))));
    if nrmref<=tol
        continue   % degenerate draw (ground truth itself ~0) -- not informative
    end

    try
        Cop = mtimes_tensopvar(Aop,Bop);
    catch ME
        error("Test_mtimes_tensopvar_zero_sweep failed: mtimes_tensopvar errored (trial=%d, seed=%d): %s", ...
            trial,70000+trial,ME.message)
    end
    % Own dummy name per factor for shared input variables (see NOTES)
    v_dum = cell(1,d); degs_j = zeros(1,N); v_full = 1;
    for j=1:d
        is_dep_j = depmat.in(j,:);
        degs_j = degs_j + is_dep_j;
        pvars_new = pvars;
        for k=find(is_dep_j)
            varname_k = pvar2varname(pvars(k));
            pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
        end
        if any(is_dep_j)
            v_dum{j} = subs(v{j},pvars(is_dep_j),pvars_new(is_dep_j));
        else
            v_dum{j} = v{j};
        end
        v_full = kron(v_full,v_dum{j});
    end
    Cv = apply_op(Cop,v_full);

    err = cleanpoly(ABv-Cv,1e-10);
    if max(max(abs(err.C)))>tol*max(1,nrmref)
        error("Test_mtimes_tensopvar_zero_sweep failed: result does not match ground truth (trial=%d, seed=%d).", ...
            trial,70000+trial)
    end
    nchecked = nchecked+1;
end

if nchecked==0
    error("Test_mtimes_tensopvar_zero_sweep: every trial was skipped or degenerate -- no trial actually exercised the check.")
end
fprintf('Test_mtimes_tensopvar_zero_sweep passed (%d/%d trials checked against ground truth, %d skipped on construction errors).\n', ...
    nchecked,ntrials,nskipped);
