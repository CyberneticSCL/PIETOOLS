function C = plus(A,B)
% C = PLUS(A,B) returns the 'tensopmat' object C representing the sum of
% the tensor-PI operators defined by the 'tensopmat' objects A and B.
%
% INPUTS
% - A:      m x n 'tensopmat' object (or 'tensopvar', 'intop', 'sopvar',
%           or 0);
% - B:      m x n 'tensopmat' object (or 'tensopvar', 'intop', 'sopvar',
%           or 0);
%
% OUTPUTS
% - C:      m x n 'tensopmat' object representing the sum of the
%           tensor-PI operators defined by A and B;
%
% NOTES
% - If A or B is 0, the other is returned.
% - A block of type 'double', 'polynomial' or 'dpvar' can be added to a
%   block of such a type, or to a 'tensopvar' block acting on a constant
%   monomial, which is converted to a 'polynomial'.
%

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

% Allow for addition with 0
if isequal(A,0)
    C = B;
    return
elseif isequal(B,0)
    C = A;
    return
end

% Check that the inputs are appropriate
if isa(A,'sopvar')
    A = tensopvar(A);
end
if isa(A,'tensopvar') || isa(A,'intop')
    A = tensopmat(A);
end
if isa(B,'sopvar')
    B = tensopvar(B);
end
if isa(B,'tensopvar') || isa(B,'intop')
    B = tensopmat(B);
end
if ~isa(A,'tensopmat') || ~isa(B,'tensopmat')
    error("Addition of 'tensopmat' objects with non-'tensopmat' objects is not supported.")
end
if ~A.chkval && ~isvalid(A)
    error("Operation not supported: The first 'tensopmat' object is not properly declared.")
else
    A.chkval = true;
end
if ~B.chkval && ~isvalid(B)
    error("Operation not supported: The second 'tensopmat' object is not properly declared.")
else
    B.chkval = true;
end

% Make sure the matrix dimensions of the operators match
if ~isequal(A.dims,B.dims)
    error("Matrix dimensions of the operators must match.")
end
% Blocks of type double/polynomial/dpvar are added to each other, or to a tensopvar that
% acts on a constant monomial, converted to a polynomial (tensopvar2polynomial). Rows
% with only such blocks map to the variables appearing in them, so express
% A, B in terms of the same variables, and merge the output spaces of such
% rows, as for a sum of polynomials (p(s)+2)
is_nat_A = native_mask(A.ops);    is_nat_B = native_mask(B.ops);
if any(is_nat_A(:)) || any(is_nat_B(:))
    for k=find(is_nat_A(:) & ~cellfun('isempty',B.ops(:)) & ~is_nat_B(:))'
        [Mk,is_ok] = tensopvar2polynomial(B.ops{k});
        if ~is_ok
            error("A block of type 'double', 'polynomial' or 'dpvar' can only be added to an operator that acts on R, a constant monomial, and can be converted to a 'polynomial', not to the operator in block %d of the second 'tensopmat'.",k)
        end
        B.ops{k} = Mk;
    end
    for k=find(is_nat_B(:) & ~cellfun('isempty',A.ops(:)) & ~is_nat_A(:))'
        [Mk,is_ok] = tensopvar2polynomial(A.ops{k});
        if ~is_ok
            error("A block of type 'double', 'polynomial' or 'dpvar' can only be added to an operator that acts on R, a constant monomial, and can be converted to a 'polynomial', not to the operator in block %d of the first 'tensopmat'.",k)
        end
        A.ops{k} = Mk;
    end
    [A,B] = common_vars(A,B);
    [A,B] = merge_out_spaces(A,B);
end
% Make sure the input and output variables match
if ~isequal(A.vars,B.vars) || ~isequal(A.dom,B.dom)
    error("Spatial variables and domains of operators must match.")
end
% Make sure the individual operators depend on the same variables
if ~isequal(A.depmat,B.depmat)
    error("Input and output function spaces of the operator in each block must match.")
end

% The operators defining the sum are given by the sum of the
% operators
C = A;
for k=1:numel(A.ops)
    if isempty(A.ops{k})
        % Assume empty operator is 0
        C.ops{k} = B.ops{k};
    elseif ~isempty(B.ops{k})
        C.ops{k} = A.ops{k} + B.ops{k};
    end
end
C.chkval = true;

end
