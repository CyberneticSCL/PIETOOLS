% Test_polyopvar_dpvar_coeffs.m -- checks 'dpvar' (decision-variable-
% carrying) content in the coefficients of a 'polyopvar':
% (1) constant2polyopvar of a dpvar constant returns the dpvar as the
%     coefficient of the degree-0 monomial;
% (2) scaling a state by a dpvar scalar (a_mult.*x, a_mult*x) is not
%     supported (the factors of a 'tensopvar' cannot depend on decision
%     variables), and gives a clear error;
% (3) a dpvar times a functional with a constant term, a_mult*(1-x'*x),
%     gives a 'dpvar' block for the constant and an 'intop' for x'*x.
% No numeric oracle exists for dpvar coefficients, so only structure and
% error messages are checked.

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 11/24/2026: Part (1) checks the dpvar constant, now supported;
% DJ, 11/24/2026: Part (2) expects the error of a dpvar times a state;
%   part (3) added;

% NOTES:
% - The random generators rand_poly/rand_sopvar never draw dpvar
%   coefficients, so no other polyopvar test covers dpvar content.

for trial = 1:15
    rng(14000+trial);

    pvar s1 s2
    x = polyopvar('x',[s1,s2],[0,1]);

    %% (1) constant2polyopvar of a genuine dpvar constant: the dpvar itself
    %% is the coefficient of the degree-0 monomial
    a_const = rand_dpvar([1,1],{'coeff1'},{},0,1);
    F = constant2polyopvar(a_const);
    if ~isequal(size(F),[1,1]) || any(F.Z(:)) || ~isa(F.C.ops{1},'dpvar')
        error("Test_polyopvar_dpvar_coeffs failed: constant2polyopvar(dpvar) is not a degree-0 monomial with a 'dpvar' coefficient (trial=%d).",trial)
    end
    D = F.C.ops{1} - a_const;
    if any(D.C(:))
        error("Test_polyopvar_dpvar_coeffs failed: the coefficient of constant2polyopvar(dpvar) differs from the dpvar (trial=%d).",trial)
    end

    %% (2) dpvar scalar times a state: a clear error, for .* and *
    a_mult = rand_dpvar([1,1],{'coeff2'},{},0,1);
    for k=1:2
        msg = '';
        try
            if k==1,    a_mult.*x;      else,   a_mult*x;   end %#ok<VUNUS>
        catch ME
            msg = ME.message;
        end
        if ~contains(msg,'is not supported')
            error("Test_polyopvar_dpvar_coeffs failed: a dpvar times a state gave '%s' instead of the error (trial=%d).",msg,trial)
        end
    end

    %% (3) dpvar times a functional with a constant: 'dpvar' block and 'intop'
    z = a_mult*(1 - x'*x);
    cls = cellfun(@class,z.C.ops,'UniformOutput',false);
    if ~isequal(sort(cls),{'dpvar','intop'})
        error("Test_polyopvar_dpvar_coeffs failed: a_mult*(1-x'*x) has blocks %s, expected a 'dpvar' and an 'intop' (trial=%d).",strjoin(cls,','),trial)
    end
end

disp('Test_polyopvar_dpvar_coeffs passed (structural checks only).')
