% Test_polyopvar_combine_terms.m -- verify combine_terms.m merges
% duplicate-monomial rows (summing their coefficients) while preserving
% the represented polynomial's value, and shrinks the monomial basis. The
% monomial row is duplicated by hand to force a nontrivial merge; oracle:
% the value of the duplicated polynomial is 2x that of the original.

% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

for trial = 1:20
    rng(4000+trial);

    pvar s1 s2
    x = polyopvar('x',[s1,s2],[0,1]);
    z = x.*x;                 % single-row Z = [1,1]
    if size(z.Z,1)~=1
        error("Test_polyopvar_combine_terms failed: setup assumption violated (trial=%d).",trial)
    end

    % duplicate the monomial row and its coefficient block
    z2 = z;
    z2.Z = [z.Z;z.Z];
    z2.C.ops = [z.C.ops, z.C.ops];
    z2.C.depmat.in = [z.C.depmat.in; z.C.depmat.in];
    z2.C.dims.in = [z.C.dims.in; z.C.dims.in];
    z2.chkval = false;
    z2.C.chkval = false;

    z3 = combine_terms(z2);
    if size(z3.Z,1)~=1
        error("Test_polyopvar_combine_terms failed: duplicate rows not merged (trial=%d).",trial)
    end
    if ~z3.chkval
        error("Test_polyopvar_combine_terms failed: chkval not set (trial=%d).",trial)
    end

    xv = rand_poly([1,1],{'s1','s2'},2,5);
    Y_orig = apply_poly(z,{xv});
    Y_dup  = apply_poly(z2,{xv});
    Y_merged = apply_poly(z3,{xv});

    err_dup = cleanpoly(Y_dup-2*Y_orig,1e-10);
    if max(max(abs(err_dup.C)))>tol
        error("Test_polyopvar_combine_terms failed: manual duplication setup did not read as 2x (trial=%d).",trial)
    end
    % merging must not change the value: z3 equals z2 (2*x^2), not the original x^2
    err_merged = cleanpoly(Y_merged-Y_dup,1e-10);
    if max(max(abs(err_merged.C)))>tol
        error("Test_polyopvar_combine_terms failed: combine_terms changed the polynomial's value (trial=%d).",trial)
    end
end

disp('Test_polyopvar_combine_terms passed.')
