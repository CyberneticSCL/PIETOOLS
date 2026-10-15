% Test_polyopvar_dpvar_coeffs.m -- checks two known gaps for 'dpvar'
% (decision-variable-carrying) content in the coefficients of a
% 'polyopvar', by asserting the failures they cause:
% (1) constant2polyopvar of a dpvar constant errors: it dispatches on
%     isa(C,'dpvar') but then calls mat2tensopvar -> mat2sopvar, which
%     reject anything but 'double'/'polynomial';
% (2) scaling a state by a dpvar scalar (a_mult.*x) succeeds and promotes
%     the leaf to an 'sdopvar', but apply_poly cannot evaluate it:
%     'sdopvar' has no apply_op method.
% No numeric oracle exists for dpvar coefficients, so only structure and
% error messages are checked.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - The test fails once either gap is closed, so that it does not stay
%   green silently; it should then be extended to check the numeric
%   value.
% - Gap (1) has its root cause in mat2sopvar.m of the read-only PIETOOLS
%   repo; gap (2) needs an apply_op counterpart for 'sdopvar'.
% - The random generators rand_poly/rand_sopvar never draw dpvar
%   coefficients, so no other polyopvar test covers dpvar content.

for trial = 1:15
    rng(14000+trial);

    pvar s1 s2
    x = polyopvar('x',[s1,s2],[0,1]);

    %% (1) constant2polyopvar of a genuine dpvar constant
    a_const = rand_dpvar([1,1],{'coeff1'},{},0,1);
    try
        constant2polyopvar(a_const);
        error("Test_polyopvar_dpvar_coeffs failed: expected the documented constant2polyopvar(dpvar) failure, but it succeeded (trial=%d) -- the maintainer may have fixed gap (1) in this test's header; update the test (and re-check whether apply_poly can now evaluate the result, see part (2)).",trial)
    catch ME
        if contains(ME.message,'Test_polyopvar_dpvar_coeffs')
            rethrow(ME);
        elseif ~contains(ME.message,"'double' or 'polynomial'")
            rethrow(ME);
        end
    end

    %% (2) dpvar scalar times a state: construction succeeds, evaluation is blocked
    a_mult = rand_dpvar([1,1],{'coeff2'},{},0,1);
    z = a_mult.*x;
    if ~isa(z,'polyopvar') || ~isequal(size(z),[1,1])
        error("Test_polyopvar_dpvar_coeffs failed: a_mult.*x did not return a 1x1 polyopvar object (trial=%d).",trial)
    end
    leaf = z.C.ops{1}.ops{1,1};
    if ~isa(leaf,'sdopvar')
        error("Test_polyopvar_dpvar_coeffs failed: expected dpvar*sopvar to promote to 'sdopvar' (got '%s', trial=%d) -- the class-level promotion behavior this test relies on may have changed; update accordingly.",class(leaf),trial)
    end

    xv = rand_poly([1,1],{'s1','s2'},2,4);
    try
        apply_poly(z,{xv});
        error("Test_polyopvar_dpvar_coeffs failed: expected the documented apply_op 'sdopvar'-dispatch failure, but evaluation succeeded (trial=%d) -- gap (2) in this test's header may have been fixed upstream; if so, extend this test to verify the numeric value (e.g. against dpvar2poly/rand_dpvar's own decision-variable convention) instead of only checking it runs.",trial)
    catch ME
        if contains(ME.message,'Test_polyopvar_dpvar_coeffs')
            rethrow(ME);
        % 'sdopvar' has no apply_op method (undefined-function error)
        elseif ~contains(ME.message,"type 'sdopvar'")
            rethrow(ME);
        end
    end
end

disp('Test_polyopvar_dpvar_coeffs passed (structural/gap check only -- see header note on dpvar-coefficient coverage being blocked upstream).')
