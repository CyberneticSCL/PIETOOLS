% Test_mtimes_tensopvar_S2_extra_var.m -- Tests the intop*tensopvar         % DJ (CC), 10/21/2026
% composition when Bop introduces an "S2-type" output variable: one that    % DJ (CC), 10/21/2026
% appears in Bop.vars.out but has NO corresponding entry in Bop.vars.in (a
% newly introduced output, as opposed to a variable shared/passed through
% between the input and output of Bop).
%
% Oracle: apply_op(Bop) followed by apply_op(Aop), independent of the
% composition.
%
% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/21/2026: Renamed from Test_mtimes_tensopvar_S1_extra_var;
%   S1 (input only) and S2 (output only) swapped, to match 'sopvar' and
%   mtimes_tensopvar, so that "S2-type" is now the output-only variable;
%
% NOTES
% - The mirror "S1-type" case (variable in Bop.vars.in only) is covered by  % DJ (CC), 10/21/2026
%   Test_mtimes_tensopvar_shared_input_var.m.
% - Deliberately small (d==1, N3<=5): the case needs no more.

tol = 1e-8;
nchecked = 0;
ntrials = 0;
for N3 = [3,4,5]
    for seed_ofst = [0,100,200]
        ntrials = ntrials+1;
        rng(50400+N3+seed_ofst);

        d = 1;
        shared_vars = arrayfun(@(k) sprintf('sShared%d',k),1:N3,'UniformOutput',false);
        vars = struct('in',{shared_vars},'out',{[shared_vars,{'sExtraOut'}]});
        depmat = struct('in',true(d,N3),'out',true(d,N3+1));

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
%           error("Test_mtimes_tensopvar_S1_extra_var failed: N3=%d case does not match ground truth (seed_ofst=%d).",N3,seed_ofst) % DJ (CC), 10/21/2026 (was)
            error("Test_mtimes_tensopvar_S2_extra_var failed: N3=%d case does not match ground truth (seed_ofst=%d).",N3,seed_ofst) % DJ (CC), 10/21/2026
        end
        nchecked = nchecked+1;
    end
end
if nchecked==0
%   error("Test_mtimes_tensopvar_S1_extra_var: every trial was skipped (all zero ground truth) -- test exercised nothing.") % DJ (CC), 10/21/2026 (was)
    error("Test_mtimes_tensopvar_S2_extra_var: every trial was skipped (all zero ground truth) -- test exercised nothing.") % DJ (CC), 10/21/2026
end
%fprintf('Test_mtimes_tensopvar_S1_extra_var passed: %d/%d trials checked against ground truth.\n',nchecked,ntrials); % DJ (CC), 10/21/2026 (was)
fprintf('Test_mtimes_tensopvar_S2_extra_var passed: %d/%d trials checked against ground truth.\n',nchecked,ntrials); % DJ (CC), 10/21/2026
