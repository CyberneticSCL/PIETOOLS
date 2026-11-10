function [logval, max_error, residual] = polyopvar_equal(F1, F2, tol)
    % [LOGVAL,MAX_ERROR,RESIDUAL] = POLYOPVAR_EQUAL(F1,F2,TOL) tests
    % whether two 'polyopvar' objects represent the same FDP.
    %
    % Direct comparison of the object fields is generally insufficient,
    % since equivalent FDPs may store repeated distributed monomials in
    % different orders. This function therefore forms F1-F2, canonicalizes
    % equivalent distributed monomials using clean_poly, and compares all
    % remaining coefficients with a specified tolerance.
    %
    % INPUTS
    % - F1, F2:     'polyopvar' objects to compare;
    % - tol:        (optional) nonnegative scalar coefficient tolerance.
    %               Defaults to 1e-10;
    %
    % OUTPUTS
    % - logval:     logical value equal to true when the maximum absolute
    %               coefficient in the canonical residual is at most tol.
    %               False if F1 and F2 have different sizes;
    % - max_error:  maximum absolute coefficient in the residual F1-F2.
    %               Inf if F1 and F2 have different sizes;
    % - residual:   canonical 'polyopvar' representation of F1-F2. Empty
    %               if F1 and F2 have different sizes;
    %
    % NOTES
    % - For dpvar-valued objects, this function tests algebraic equality as
    %   affine expressions in the decision variables. It does not use any
    %   equality constraints stored separately in a PIESOS program.
    % - To compare objects after solving a PIESOS program, first substitute
    %   the solution using piesos_getsol and then call this function.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PIETOOLS - polyopvar_equal
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
    % If you modify this code, document all changes carefully and
    % include date authorship, and a brief description of modifications.
    %
    % CR, 10/07/2026: Initial coding.
    % DJ, 10/15/2026: Updated to use new polyopvar functions;
    % DJ, 11/09/2026: Sum the terms (rows of ops) of a 'tensopvar' 
    %                   coefficient before testing for zero;

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
    residual = clean_poly(F1-F2);                                           % DJ, 10/14/2026

    max_error = 0;
    for term_idx = 1:numel(residual.C.ops)
        coefficient = residual.C.ops{term_idx};

        term_error = max_abs_coef(coefficient);
        max_error  = max(max_error,term_error);
    end

    logval = max_error <= tol;

end


function m = max_abs_coef(coefficient)                                      % DJ, 10/14/2026
    % m = max_abs_coef(coefficient) returns the maximal absolute value of the
    % coefficients (of the kernels, for an 'intop', of the parameters, for a
    % 'sopvar') of the operator or polynomial 'coefficient'.
    %
    % INPUT
    % - coefficient:    'intop', 'tensopmat', 'tensopvar', 'sopvar', 'dpvar',
    %                   'polynomial', numeric array, or cell of these;
    %
    % OUTPUT
    % - m:              nonnegative scalar, the maximal absolute coefficient;

    m = 0;
    if isempty(coefficient)
        return
    elseif isa(coefficient,'intop')
        % The params of an 'intop' are coefficients in a unified monomial
        % basis, not kernels: use the kernels
        m = max_abs_coef(kernels(coefficient));
    elseif isa(coefficient,'tensopvar') && ...
            all(cellfun(@(c) isa(c,'sopvar'),coefficient.ops(:)))
        % Terms of a sum of tensor products must be added, not compared     % DJ, 11/09/2026
        m = max_abs_terms(coefficient.ops);
    elseif isa(coefficient,'tensopmat') || isa(coefficient,'tensopvar')
        m = max_abs_coef(coefficient.ops);
    elseif isa(coefficient,'sopvar')
        m = max_abs_coef(coefficient.params);
    elseif iscell(coefficient)
        for k=1:numel(coefficient)
            m = max(m,max_abs_coef(coefficient{k}));
        end
    elseif isa(coefficient,'dpvar')
        m = max_abs_coef(coefficient.C);
    elseif isa(coefficient,'polynomial')
        m = max_abs_coef(coefficient.coefficient);
    elseif isnumeric(coefficient)
        m = max([0;abs(full(coefficient(:)))]);
    else
        error('polyopvar_equal:UnsupportedCoefficient', ...
            ['Cannot compare a polyopvar coefficient of class ', ...
             '"%s".'],class(coefficient));
    end

end



function m = max_abs_terms(ops)                                             % DJ), 11/09/2026
    % m = max_abs_terms(ops) returns the maximal absolute coefficient of the
    % tensor T = sum_i ops{i,1} (x) ... (x) ops{i,d}, with the terms i merged
    % wherever possible so that cancelling terms give 0.
    %
    % Two terms agreeing up to a scalar in all but one factor j are replaced
    % by one term, with factor j the (scaled) sum of their factors j. This
    % is exact for the tensor T of the parameters, whose maximal entry is
    % the product of the maximal entries of the factors of the last
    % remaining term; but it is only sufficient for equality of the
    % operators when the product is pointwise (otype 1), where different T
    % can give the same operator, and it does not detect a sum that
    % vanishes only through three or more dependent terms.
    %
    % INPUT
    % - ops:    K x d cell of 'sopvar' objects, the d factors of K terms;
    %
    % OUTPUT
    % - m:      nonnegative scalar, the maximal absolute coefficient of T;

    ztol = 1e-12;                       % as in clean_poly
    [K,d] = size(ops);

    % Split each factor into a unit-max factor Rn and its scale s, so that
    % proportional factors compare equal (up to sign)
    s = zeros(K,d);
    Rn = ops;
    for i=1:K
        for j=1:d
            s(i,j) = max_abs_coef(ops{i,j});
            if s(i,j)>ztol
                Rn{i,j} = (1/s(i,j))*ops{i,j};
            end
        end
    end
    % A term with a zero factor is zero
    alive = all(s>ztol,2);
    w = prod(s,2);                      % term i = w(i) * Rn{i,1} (x) ... (x) Rn{i,d}

    % Merge term k into term i if they agree in all factors but (at most) one
    changed = true;
    while changed
        changed = false;
        for i=1:K
            for k=i+1:K
                if ~alive(i) || ~alive(k)
                    continue
                end
                % Rn{k,m} = sg(m)*Rn{i,m}, or sg(m)=0 if not proportional
                sg = zeros(1,d);
                for m=1:d
                    if max_abs_coef(Rn{i,m}-Rn{k,m})<=ztol
                        sg(m) = 1;
                    elseif max_abs_coef(Rn{i,m}+Rn{k,m})<=ztol
                        sg(m) = -1;
                    end
                end
                j = find(sg==0);
                if numel(j)>1
                    continue
                elseif isempty(j)
                    j = 1;
                end
                % Add the factors j, scaling by the other factors of k
                sg(j) = 1;
                Q = w(i)*Rn{i,j} + (w(k)*prod(sg))*Rn{k,j};
                w(i) = max_abs_coef(Q);
                alive(k) = false;
                changed = true;
                if w(i)>ztol
                    Rn{i,j} = (1/w(i))*Q;
                else
                    alive(i) = false;   % terms cancelled
                end
            end
        end
    end

    % Factors Rn have unit max, so the max entry of term i is w(i)
    m = max([0;w(alive)]);

end
