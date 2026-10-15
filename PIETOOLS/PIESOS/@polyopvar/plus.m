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
% - Addition of polynomials depending on different spatial variables is
%   not supported: the output dimensions must match, e.g. p1(x)(s1) +
%   p2(y)(s2) does not work, as we would have to augment p1(x)(s1,s2) and
%   p2(y)(s1,s2), which is not implemented.

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