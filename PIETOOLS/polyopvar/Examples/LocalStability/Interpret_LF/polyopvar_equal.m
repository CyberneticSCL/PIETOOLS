function [logval, max_error, residual] = polyopvar_equal(F1, F2, tol)
    % [logval,max_error,residual] = polyopvar_equal(F1,F2,tol) tests
    % whether two polyopvar objects represent the same FDP.
    %
    % Direct comparison of the object fields is generally insufficient,
    % since equivalent FDPs may store repeated distributed monomials in
    % different orders. This function therefore forms F1-F2, canonicalizes
    % equivalent distributed monomials using piesos_combine_terms, and
    % compares all remaining coefficients with a specified tolerance.
    %
    % INPUTS
    % - F1, F2      polyopvar objects to compare.
    % - tol         Optional nonnegative scalar coefficient tolerance.
    %               The default value is 1e-10.
    %
    % OUTPUTS
    % - logval      Logical value equal to true when the maximum absolute
    %               coefficient in the canonical residual is at most tol.
    % - max_error   Maximum absolute coefficient in the residual F1-F2.
    % - residual    Canonical polyopvar representation of F1-F2.
    %
    % NOTES
    % - For dpvar-valued objects, this function tests algebraic equality as
    %   affine expressions in the decision variables. It does not use any
    %   equality constraints stored separately in a PIESOS program.
    % - To compare objects after solving a PIESOS program, first substitute
    %   the solution using piesos_getsol and then call this function.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - polyopvar_equal.m
    %
    % Copyright (C) 2026 PIETOOLS Team
    %
    % This program is free software; you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation; either version 2 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program; if not, write to the Free Software
    % Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % If you modify this code, document all changes carefully and include date
    % authorship, and a brief description of modifications.
    %
    % CR, 10/07/2026: Initial coding.

    narginchk(2,3);

    if nargin < 3
        tol = 1e-10;
    end

    if ~isa(F1,'polyopvar') || ~isa(F2,'polyopvar')
        error('polyopvar_equal:InvalidInput', ...
            'F1 and F2 must both be polyopvar objects.');
    end

    if ~isnumeric(tol) || ~isscalar(tol) || ~isreal(tol) || ...
            ~isfinite(tol) || tol < 0
        error('polyopvar_equal:InvalidTolerance', ...
            'tol must be a finite, nonnegative, real scalar.');
    end

    if ~isequal(size(F1),size(F2))
        logval   = false;
        max_error = Inf;
        residual = [];
        return
    end

    % Canonicalization accounts for equivalent distributed monomials such
    % as x(t1)*x(t2) and x(t2)*x(t1).
    residual = piesos_combine_terms(F1-F2);

    max_error = 0;
    for term_idx = 1:numel(residual.C.ops)
        coefficient = residual.C.ops{term_idx};

        if isempty(coefficient)
            continue
        elseif isa(coefficient,'intop')
            coefficient = coefficient.params;
        end

        if isa(coefficient,'dpvar')
            values = coefficient.C;
        elseif isa(coefficient,'polynomial')
            values = coefficient.coefficient;
        elseif isnumeric(coefficient)
            values = coefficient;
        else
            error('polyopvar_equal:UnsupportedCoefficient', ...
                ['Cannot compare a polyopvar coefficient of class ', ...
                 '"%s".'],class(coefficient));
        end

        if ~isempty(values)
            term_error = max(abs(full(values(:))));
            max_error  = max(max_error,term_error);
        end
    end

    logval = max_error <= tol;

end
