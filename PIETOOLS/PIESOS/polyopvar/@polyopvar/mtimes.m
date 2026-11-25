function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'polyopvar' object C representing the product
% A*B of two distributed polynomials, or of a constant or an operator with
% a distributed polynomial.
%
% INPUTS
% - A:     'polyopvar' object; or a 'double' matrix, a scalar 'polynomial'
%          (a constant, or a function of the output variables of B) or
%          'dpvar', or an operator ('sopvar', 'tensopvar' or
%          'tensopmat') acting on the output of B, such that A*B.C is
%          defined;
% - B:     'polyopvar' object, or a scalar 'double', 'polynomial' or
%          'dpvar';
%
% OUTPUTS
% - C:     'polyopvar' object representing A*B;
%
% NOTES
% - The product of two 'polyopvar' objects is only supported if one of
%   them is scalar, and is then the elementwise product A.*B. The product
%   of a functional with a function of the spatial variables, e.g.
%   x*(x'*x), is not supported.
% - A decision variable ('dpvar') c may multiply functionals and the
%   constant term, e.g. c*(1-x'*x), but not a state, e.g. c*x, also
%   without spatial variables; use c as a weight instead, x'*c*x.
% - A 'polyopvar' object A can only be right-multiplied by a scalar B,
%   and then A*B = B*A.
% - A scalar function p(s) of the output variables multiplies pointwise,
%   (p*B)(x)(s) = p(s)*B(x)(s).
% - For a non-scalar 'double' matrix A, each row of A*B is a linear
%   combination of the rows of B. Rows of B with different output spaces
%   (e.g. R and L2[s]) are first lifted to a common one, e.g. an element
%   of R to a function of s that is constant in s. A monomial acted on by
%   an operator with several factors (e.g. x0*u) is not supported.
% - If A is an adjoint (A=U'*W), A*B is
%   - for a 'polyopvar' B, the inner product innerprod(U,W,B), a scalar
%     'polyopvar', so that U'*B = innerprod(U,B);
%   - for a 'double', 'polynomial' or 'dpvar' B with as many rows as A has
%     columns, or a scalar, the adjoint U'*W*B, with weight W*B, so that
%     U'*P*V = innerprod(U,P,V). A scalar weight P represents P*eye(n).
%   An adjoint as the right factor, e.g. U*V', is only supported for a
%   constant scalar A, with A*U'*W = U'*(W*A), e.g. 2*u'*f.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
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
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

% Make sure the inner dimensions match
if size(A,2)~=size(B,1) && ~all(size(A)==1) && ~all(size(B)==1)
    error("Incorrect dimensions for matrix multiplication.")
end

% Handle an adjoint separately: it is not a polynomial in the states
if isa(B,'polyopvar') && B.is_adjoint
    if ~isa(A,'polyopvar') && all(size(A)==1) && (isnumeric(A) || ((isa(A,'polynomial') || isa(A,'dpvar')) && isempty(A.varname)))
        % A constant scalar commutes with the adjoint: A*U'*W = U'*(W*A)
        C = B*A;
        return
    end
    error("Multiplication with an adjoint as the right factor, e.g. u*v', is not supported.")
elseif isa(A,'polyopvar') && A.is_adjoint
    % The polynomial U of A=U'*W, without the weight W
    Au = A;     Au.is_adjoint = false;      Au.weight = 1;
    if isa(B,'polyopvar')
        % <U,W*B>: innerprod checks the dimensions of W, as the check above
        % skips a scalar A
        C = innerprod(Au,A.weight,B);
    elseif (isnumeric(B) || isa(B,'polynomial') || isa(B,'dpvar')) && (all(size(B)==1) || size(B,1)==size(A,2))
        % U'*W*B stores the weight W*B, for the inner product with a later
        % factor. Check the dimensions of B here, as the check above skips
        % a scalar A
        C = A;
        C.weight = A.weight*B;
    else
        error("An adjoint can only be multiplied with a 'polyopvar' object, or with a scalar or matrix of type 'double', 'polynomial' or 'dpvar' with as many rows as the adjoint has columns.")
    end
    return
end

% Distinguish three cases
if isa(A,'polyopvar') && isa(B,'polyopvar')
    % Support multiplication of polynomials only if one is scalar for now
    if all(size(A)==1) || all(size(B)==1)
        C = A.*B;
        return
    else
        error("Multiplication of matrix-valued distributed polynomials is currently not supported.")
    end
elseif ~isa(A,'polyopvar')
    if isa(A,'polynomial') && all(size(A)==1) && ~isempty(A.varname)
        % A scalar function of the spatial variables of B multiplies B
        % pointwise, as a degree-0 polynomial on the domain of B
        C = constant2polyopvar(A,B.vars.out,B.dom.out).*B;
        return
    end
    % A column times a scalar is the elementwise product, as for matrices
    % (the operators of B may be multi-factor, which '*' does not support)
    if isnumeric(A) && size(A,1)>1 && size(A,2)==1 && all(size(B)==1)
        C = A.*B;
        return
    end
    % A decision variable multiplies the constant term and functionals
    % ('intop' blocks), not a state, also without spatial variables
    if isa(A,'dpvar') && any(cellfun(@(op) ~isempty(op) && ~isa(op,'intop'),B.C.ops(:,full(any(B.Z,2)))),'all')
        error("Multiplication of a decision variable ('dpvar') with a distributed polynomial in the states, e.g. c*x, is not supported; multiply a functional instead, e.g. c*(x'*x), or use c as a weight, e.g. x'*c*x.")
    end
    % Multiplication of constant with polynomial
    % A row of A*B combines the rows of B, so these must map to the same
    % output space: lift them to a common one if they differ.
    % An operator A acts on the output space of each row as it is: no lift.
    if isa(A,'double') && ~all(size(A)==1)
        if ~B.chkval
            B = initialize(B);
        end
        if any(any(diff(B.depmat.out,1,1)))
            B = lift_output(B,B.vars,B.dom,[],max(B.depmat.out,[],1));
        end
    end
    C = B;
    % The output space of C is set by the assignment of C.C
    C.C = A*B.C;
    return
elseif ~isa(B,'polyopvar')
    % Multiplication of polynomial with constant
    if isnumeric(B) && all(size(B)==1)
        C = A;
        C.C = A.C*B;
        return
    elseif (isa(B,'polynomial') || isa(B,'dpvar')) && all(size(B)==1)
        % A scalar coefficient commutes with the polynomial
        C = B*A;
        return
    else
        error("A 'polyopvar' object can only be right-multiplied by a scalar.")
    end
end

end