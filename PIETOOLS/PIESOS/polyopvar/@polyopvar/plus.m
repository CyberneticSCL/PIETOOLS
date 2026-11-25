function C = plus(A,B)
% C = PLUS(A,B) returns the 'polyopvar' object C representing the sum A+B
% of two distributed polynomials, or of a distributed polynomial and a
% constant.
%
% INPUTS
% - A:      m x n 'polyopvar' object, or a 'double' or 'dpvar' column
%           vector (a constant), or a 'polynomial' column vector in the
%           output variables of B (a known function);
% - B:      m x n 'polyopvar' object, or a constant or known function, as
%           A;
%
% OUTPUTS
% - C:      m x n 'polyopvar' object representing A+B;
%
% NOTES
% - Adding the scalar 0 returns the other input unchanged.
% - A 'dpvar' can only be added to a polynomial without output
%   variables, e.g. V - c for a functional V.
% - Polynomials with different output spaces are lifted to a common one,
%   along whose new variables they are constant. E.g. the sum of
%   p1(x)(s1) and p2(y)(s2) is the function p1(x)(s1) + p2(y)(s2) of
%   (s1,s2), and the sum of a polynomial in an ODE state (output in R)
%   and one in a PDE state is a function of the spatial variable.
% - A scalar-valued polynomial (or constant) is added to each row of a
%   vector-valued one. Otherwise, the sizes must match.
% - An adjoint U1'*W can only be added to an adjoint U2'*W of the same
%   size, with the same weight W, giving (U1+U2)'*W.

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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% DJ, 11/24/2026: Initial coding;

% Allow addition of polynomial with 0
if isnumeric(A) && isequal(A,0)
    C = B;
    return
elseif isnumeric(B) && isequal(B,0)
    C = A;
    return
end

% A decision variable is only added as constant term of a polynomial
% without output variables (e.g. a functional)
if (isa(A,'dpvar') && ~isempty(B.vars.out)) || (isa(B,'dpvar') && ~isempty(A.vars.out))
    error("A 'dpvar' can only be added to a distributed polynomial without output variables, such as a functional u'*u.")
end
% Convert constant values to 'polyopvar' class; a function of spatial
% variables is defined on the domain of the other term
if isa(A,'double') || isa(A,'polynomial') || isa(A,'dpvar')
    A = constant2polyopvar(A,B.vars.out,B.dom.out);
elseif ~isa(A,'polyopvar')
    error("Addition of 'polyopvar' objects with non-'polyopvar' objects is not supported.")
end
if isa(B,'double') || isa(B,'polynomial') || isa(B,'dpvar')
    B = constant2polyopvar(B,A.vars.out,A.dom.out);
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
    elseif ~isequal(A.vars.out,B.vars.out) || ~isequal(A.dom.out,B.dom.out) || ~isequal(A.depmat.out,B.depmat.out)
        % (U1+U2)'*W lifts U1 and U2 to a common output space, which would
        % scale their inner products by the length of the added domains
        error("Adjoints can only be added if they have the same output space; add the inner products instead, e.g. x'*v + y'*v.")
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
[A,B] = common_outputs(A,B);

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
