% Test_mtimes_tensopvar_S3_mixed_deps.m -- Tests the intop*tensopvar
% composition when Bop has N3=3 shared (pass-through) variables, but its
% factors do NOT all depend on the same ones, e.g. factor 1 on (sA,sB) and
% factor 2 on (sA,sB,sC).
%
% Oracle: apply_op(Bop) followed by apply_op(Aop) on concrete polynomials,
% independent of the composition.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES
% - Tests the handling of N3-D parameter cells of different ndims in
%   misc/otimes_params.m>otimes_quadPolyCell: MATLAB drops trailing
%   singleton dimensions of the N3-D parameter cells, so factors that do
%   not depend on the last shared variable differ in ndims. Needs N3>=3,
%   different shared-variable sets across the factors, and d>=2 (otimes is
%   only called to combine factors).
% - Each factor's copy of a shared input variable gets a distinct dummy
%   name (_dum<j>) before the copies are combined into the input of Cop.

tol = 1e-8;
vars = struct('in',{{'sA','sB','sC'}},'out',{{'sA','sB','sC'}});
dom = [0,1];
dims = [2,1; 1,2];
N = numel(vars.in);     d = size(dims,1);

% Each case lists, per factor (rows), on which of (sA,sB,sC) it depends
%   in:  input variables of the factor;     out: output variables
cases = struct('name',{},'in',{},'out',{});
cases(end+1) = struct('name','x1(sA,sB) & x2(sA,sB,sC)', ...
    'in',logical([1,1,0; 1,1,1]),'out',logical([1,1,0; 1,1,1]));
cases(end+1) = struct('name','x1(sA,sB,sC) & x2(sA,sB)', ...
    'in',logical([1,1,1; 1,1,0]),'out',logical([1,1,1; 1,1,0]));
cases(end+1) = struct('name','x1(sA,sB) & x2(sC)', ...
    'in',logical([1,1,0; 0,0,1]),'out',logical([1,1,0; 0,0,1]));
cases(end+1) = struct('name','x1(sC) & x2(sA,sB)', ...
    'in',logical([0,0,1; 1,1,0]),'out',logical([0,0,1; 1,1,0]));
cases(end+1) = struct('name','common output space, inputs (sA,sB) & (sA,sB,sC)', ...
    'in',logical([1,1,0; 1,1,1]),'out',true(2,3));

nchecked = 0;
for c = 1:numel(cases)
    ok = false;
    for seed = 1:5      % retry if the random draw makes the reference vanish
        rng(140000+10*c+seed);
        depmat = struct('in',cases(c).in,'out',cases(c).out);
        Bop = rand_tensopvar(dims,vars,dom,depmat,1,[2,2],1,0.5);
        pdegs = ones(1,numel(Bop.vars.out));
        Aop = rand_intop([1,Bop.matdim(1)],Bop.vars.out,Bop.dom.out,Bop.deps.out,pdegs);

        % Random input to each factor, and the same input with each factor's
        % own dummy copy of every shared variable
        pvars = polynomial(Bop.vars.in(:));
        v = cell(1,d);      v_dum = cell(1,d);
        degs_j = zeros(1,N);        v_full = 1;
        for j=1:d
            is_dep_j = cases(c).in(j,:);
            n_j = size(Bop.ops{1,j},2);
            v{j} = rand_poly([n_j,1],pvars(is_dep_j),2,50);
            degs_j = degs_j + is_dep_j;
            pvars_new = pvars;
            for k=find(is_dep_j)
                varname_k = pvar2varname(pvars(k));
                pvars_new(k) = polynomial({[varname_k{1},'_dum',num2str(degs_j(k))]});
            end
            v_dum{j} = subs(v{j},pvars(is_dep_j),pvars_new(is_dep_j));
            v_full = kron(v_full,v_dum{j});
        end

        % Ground truth: apply the operators one after the other
        Bv = apply_op(Bop,v);
        ABv = apply_op(Aop,Bv);
        nrmref = full(max(max(abs(double(ABv)))));
        if nrmref<=tol
            continue
        end

        % Composition
        Cop = mtimes_tensopvar(Aop,Bop);
        Cv = apply_op(Cop,v_full);
        err = cleanpoly(ABv-Cv,1e-10);
        if max(max(abs(err.C)))>tol*max(1,nrmref)
            error("Test_mtimes_tensopvar_S3_mixed_deps failed: case '%s' does not match ground truth (seed=%d).",cases(c).name,seed)
        end
        ok = true;
        break
    end
    if ~ok
        error("Test_mtimes_tensopvar_S3_mixed_deps: case '%s' -- every random draw gave a (near) zero reference, nothing exercised.",cases(c).name)
    end
    nchecked = nchecked+1;
end
fprintf('Test_mtimes_tensopvar_S3_mixed_deps passed (%d/%d cases checked against ground truth).\n',nchecked,numel(cases));
