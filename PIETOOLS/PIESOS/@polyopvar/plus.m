function C = plus(A,B)
% C = PLUS(A,B) returns the 'polyopvar' object C representing the sum
% of the distributed polynomials defined by the 'polyopvar' objects A and
% B.
%
% INPUTS
% - A:      m x n 'polyopvar' object, or a 'double', 'polynomial' or
%           'dpvar' array representing a constant;
% - B:      m x n 'polyopvar' object, or a 'double', 'polynomial' or
%           'dpvar' array representing a constant;
%
% OUTPUTS
% - C:      m x n 'polyopvar' object representing the sum of the
%           distributed polynomials defined by A and B;
%
% NOTES
% - Adding the scalar 0 to a 'polyopvar' object returns that object.
% - If the output spaces of the polynomials differ, they are lifted to the
%   same space by 'lift_output': the output of a polynomial is constant
%   along the output variables it is not defined in. For example, for
%   p1(x)(s1) and p2(y)(s2), the sum is the function of (s1,s2) defined
%   by p1(x)(s1) + p2(y)(s2), and the sum of a polynomial in an ODE state
%   (output in R) and one in a PDE state is a function of the spatial
%   variable.
% - A scalar-valued polynomial (or constant) can be added to a
%   vector-valued polynomial, and is then added to each of its rows.
%   Otherwise, the vector dimensions must match.
% - An adjoint (A=U1', see 'ctranspose') can only be added to an adjoint of
%   the same size with the same weight W, as then U1'*W + U2'*W =
%   (U1+U2)'*W. The adjoints of scalar polynomials are not added to
%   adjoints of vectors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus
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
% DJ, 10/15/2026: Initial coding;
% DJ, 11/03/2026: Add support for addition of polynomials with different
%                   output spaces (e.g. p1(x1)(s1) + p2(x2)(s2));
% DJ, 11/07/2026: Add adjoints with the same weight, see 'ctranspose';

% Allow addition of polynomial with 0
if isnumeric(A) && isequal(A,0)
    C = B;
    return
elseif isnumeric(B) && isequal(B,0)
    C = A;
    return
end

% Convert constant values to 'polyopvar' class
if isa(A,'double') || isa(A,'polynomial') || isa(A,'dpvar')
    A = constant2polyopvar(A);
elseif ~isa(A,'polyopvar')
    error("Addition of 'polyopvar' objects with non-'polyopvar' objects is not supported.")
end
if isa(B,'double') || isa(B,'polynomial') || isa(B,'dpvar')
    B = constant2polyopvar(B);
elseif ~isa(B,'polyopvar')
    error("Addition of 'polyopvar' objects with non-'polyopvar' objects is not supported.")
end

% Add adjoints with equal weight W: U1'*W + U2'*W = (U1+U2)'*W, by adding
% the polynomials that they are the adjoint of
if A.is_adjoint || B.is_adjoint
    % The adjoints must have the same size and weight. isequal is not
    % defined for dpvar, but A.weight-B.weight is then a dpvar, zero if
    % it has no nonzero coefficients. isequal of 'polynomial' objects is
    % elementwise
    is_same = A.is_adjoint && B.is_adjoint && all(size(A)==size(B)) && isequal(size(A.weight),size(B.weight));
    if is_same && (isa(A.weight,'dpvar') || isa(B.weight,'dpvar'))
        D = A.weight - B.weight;
        is_same = ~any(any(D.C));
    elseif is_same
        is_same = all(all(isequal(A.weight,B.weight)));
    end
    if ~is_same
        error("An adjoint can only be added to an adjoint of the same size, with the same weight.")
    end
    W = A.weight;
    A.is_adjoint = false;       A.weight = 1;
    B.is_adjoint = false;       B.weight = 1;
    C = plus(A,B);
    C.is_adjoint = true;        C.weight = W;
    return
end

% Make sure the two distributed polynomials have the same output space, so
% that the operators defining them can be added
[A,B] = common_outputs(A,B);                                                % DJ, 11/03/2026

% Make sure the sizes of the polynomials match
if any(size(A)~=size(B))
    error("Dimensions of the polynomials to add must match.")
end

% Make sure the two distributed polynomials are expressed in terms of the
% same basis of distributed monomials.
[A,B] = common_basis(A,B);

% Then, the coefficients defining the sum are given by the sum of the
% coefficients
C = A;
C.C = A.C + B.C;
C.chkval = true;
C = combine_terms(C);

end
