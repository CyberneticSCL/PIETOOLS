% Test_mtimes_tensopvar_dvar_lists.m -- Tests the intop*tensopvar
% composition when the kernels of the functional depend on DIFFERENT
% decision variables in different regions of integration, as for the
% functionals declared by piesos_quadvar. The composition adds the
% contributions of different regions (plus_spstruct in
% mtimes_tensopvar.m).
%
% Oracle: the decision variables are evaluated at random values d0 BEFORE
% composing, and the functional applied to the output of Bop (apply_op, no
% composition). The composition of the symbolic functional, evaluated at
% d0, must agree.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - All kernels of an 'intop' share one list Zd: "different decision
%   variables" means coefficients on different subsets of Zd (thin_dvars
%   zeroes the other rows, used_dvars reads the subsets).

tol = 1e-8;
ntested = 0;
q = 12;

dims = [2,1; 1,2];
vars = struct('in',{{'sShared'}},'out',{{'sOut1','sShared'}});
dom = [0,1];
depmat_cases = {logical([1,1;1,0]), logical([1,1;0,1]), logical([1,1;1,1])};

for c = 1:numel(depmat_cases)
    for nterms = 1:2
        depmat = struct('in',true(2,1),'out',depmat_cases{c});
        % The test is only meaningful if the functional has several regions, with kernels that
        % depend on different decision variables: redraw until it does (the random integrals
        % may reduce to a single region)
        for attempt = 1:50
            rng(5000+10*c+nterms+1000*attempt);
            Bop = rand_tensopvar(dims,vars,dom,depmat,nterms,[2,2],1,0.5);
            pdegs = ones(1,numel(Bop.vars.out));
            Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs,q);
            Aop = thin_dvars(Aop);
            lists = cellfun(@(P) used_dvars(P,numel(Aop.Zd)),Aop.params(:),'UniformOutput',false);
            if numel(unique(lists))>=2
                break
            end
        end
        if numel(unique(lists))<2
            error("Test_mtimes_tensopvar_dvar_lists: all kernels depend on the same decision variables (case=%d, nterms=%d).",c,nterms)
        end
        d0map = containers.Map(arrayfun(@(i) ['coeff',num2str(i)],1:q,'UniformOutput',false),num2cell(randn(q,1)));
        Aop_num = eval_intop_dpvar(Aop,d0map);

        pvars = polynomial(Bop.vars.in(:));
        v = {rand_poly([1,1],pvars,2,10), rand_poly([2,1],pvars,2,10)};
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop_num,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue
        end

        Cop = mtimes_tensopvar(Aop,Bop);
        % Both factors depend on the SAME sole input variable 'sShared': give each its
        % own dummy copy before combining (see Test_mtimes_tensopvar_shared_input_var.m)
        d = 2;  degs_j = 0;  v_full = 1;  v_dum = cell(1,d);
        for j=1:d
            degs_j = degs_j + 1;
            varname_1 = pvar2varname(pvars(1));
            pvar_new = polynomial({[varname_1{1},'_dum',num2str(degs_j)]});
            v_dum{j} = subs(v{j},pvars(1),pvar_new);
            v_full = kron(v_full,v_dum{j});
        end
        Cv = eval_dpvar_at(apply_op(Cop,v_full),d0map);
        err = cleanpoly(ABv-Cv,1e-10);
        ntested = ntested+1;
        if ~isdouble(err) || max(max(abs(double(err))))>tol*max(1,nrmref)
            error("Test_mtimes_tensopvar_dvar_lists failed: result does not match ground truth (case=%d, nterms=%d).",c,nterms)
        end
    end
end

if ntested==0
    error("Test_mtimes_tensopvar_dvar_lists: every case was skipped or degenerate -- no case actually exercised the check.")
end
fprintf('Test_mtimes_tensopvar_dvar_lists passed (%d/%d cases checked).\n',ntested,numel(depmat_cases)*2);


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function Aop = thin_dvars(Aop)
% Retain a different random subset of the decision variables in each parameter of Aop
% (the rows of the others, in the blocks of rows of each row of the matrix, are set to zero)
nd = numel(Aop.Zd);
m = size(Aop,1);
for k=1:numel(Aop.params)
    P = Aop.params{k};
    if nd==0,       continue,       end
    J = sort(randperm(nd,ceil(nd/2)));
    rows = reshape(((0:m-1)'*(nd+1) + [1,1+J]).',[],1);
    P(setdiff(1:size(P,1),rows),:) = 0;
    Aop.params{k} = P;
end
end

function s = used_dvars(P,nd)
% String with the decision variables (rows of P) that have nonzero coefficients
[i,~] = find(P);
d = mod(i-1,nd+1);          % 0 for the constant term
s = char(strjoin(string(unique(d(d>0))'),','));
end
