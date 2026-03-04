function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'tensopvar' object C representing the product
% of the multiplier operator defined by A with the tensor-PI operator
% defined by B.
%
% INPUTS
% - A:      m x p 'double', 'polynomial', or 'dpvar' object representing 
%           some multiplier operator to multiply with B;
% - B:      p x n 'tensopvar' object acting on some degree-d distributed 
%           monomial;
%
% OUTPUTS
% - C:      m x n 'tensopvar' object representing the product of the
%           multiplier A with the tensor-PI operator defined by B;
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
% DJ, 03/03/2026: Initial coding

% Check that multiplication is supported
if isa(A,'tensopvar') && isa(B,'tensopvar')
    error("Multiplication of tensor-PI operators is not supported.")
elseif isa(A,'tensopvar') && all(size(B)==1)
    % Allow right-multiplication with a scalar
    C = B*A;
    return
elseif isa(A,'tensopvar')
    error("Left-multiplication with 'tensopvar' objects is not supported.")
elseif ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'tensopvar' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end

% Check that the inner dimensions of the objects match
Adim = size(A);
Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Prohibit matrix-multiplication with actual tensor products of PI
% operators for now, since this is not quite trivial
%   (note that  A*(B1*x)(s)*(B2*x)(s) ~= (A*B1*x)(s)*(B2*x)(s))
if ~all(Adim==1) && any(C1.degmat>1)
    error("Matrix-multiplication with tensor-products of PI operators is currently not supported. Use '.*' for elementwise multiplication instead.")
end

% Multiply each of the operators stored in B by A
C = B;
for j=1:numel(C.ops)
    Cj = C.ops{j};
    if isa(Cj,'cell')
        % Multiply the first factor of each term with A
        for i=1:size(Cj,1)
            Cj{i,1} = A*Cj{i,1};
        end
    else
        % Multiply the operator acting on the jth monomial with A
        Cj = A*Cj;
    end
    C.ops{j} = Cj;
end

end