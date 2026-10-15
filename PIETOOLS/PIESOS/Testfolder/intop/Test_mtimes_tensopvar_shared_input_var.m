% Test_mtimes_tensopvar_shared_input_var.m -- Tests the intop*tensopvar
% composition when >=2 of the factors of Bop depend on the SAME named input
% variable (so their copies of it must be given distinct dummy names before
% being combined into a single test input, as in
% Test_intop_mtimes_tensopvar.m).
%
% Cases: symmetric and asymmetric per-factor dependence on a shared output
% variable, d=2, nterms in {1,2}. Oracle: apply_op(Bop) followed by
% apply_op(Aop).
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - The test input MUST give each factor its own dummy copy of the shared
%   variable (_dum<j>): kron(v{1},v{2}) with both in the one variable name
%   conflates two independent copies and mismatches the (correct)
%   composition.

tol = 1e-8;
ntested = 0;

dims = [2,1; 1,2];
vars = struct('in',{{'sShared'}},'out',{{'sOut1','sShared'}});
dom = [0,1];
depmat_cases = {logical([1,1;1,0]), logical([1,1;0,1]), logical([1,1;1,1])};

for c = 1:numel(depmat_cases)
    for nterms = 1:2
        rng(999);
        depmat = struct('in',true(2,1),'out',depmat_cases{c});
        Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],1,0.5);
        pdegs = ones(1,numel(Bop.vars.out));
        Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

        pvars = polynomial(Bop.vars.in(:));
        v = {rand_poly([1,1],pvars,2,10), rand_poly([2,1],pvars,2,10)};
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue
        end

        Cop = mtimes_tensopvar(Aop,Bop);
        % Both factors depend on the SAME sole input variable 'sShared' --
        % give each its own dummy copy before combining (see NOTES).
        d = 2; N = 1; degs_j = 0; v_full = 1; v_dum = cell(1,d);
        for j=1:d
            degs_j = degs_j + 1;
            varname_1 = pvar2varname(pvars(1));
            pvar_new = polynomial({[varname_1{1},'_dum',num2str(degs_j)]});
            v_dum{j} = subs(v{j},pvars(1),pvar_new);
            v_full = kron(v_full,v_dum{j});
        end
        Cv = apply_op(Cop,v_full);
        err = cleanpoly(ABv-Cv,1e-10);
        ntested = ntested+1;
        if max(max(abs(err.C)))>tol*max(1,nrmref)
            error("Test_mtimes_tensopvar_shared_input_var failed: result does not match ground truth (case=%d, nterms=%d).",c,nterms)
        end
    end
end

if ntested==0
    error("Test_mtimes_tensopvar_shared_input_var: every case was skipped or degenerate -- no case actually exercised the check.")
end
fprintf('Test_mtimes_tensopvar_shared_input_var passed (%d/%d cases checked).\n',ntested,numel(depmat_cases)*2);
