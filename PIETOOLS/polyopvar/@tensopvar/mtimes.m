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
% DJ, 04/22/2026: Initial coding

% Check that multiplication is supported
if ~isa(B,'tensopvar') && all(size(B)==1)
    % Allow right-multiplication with a scalar
    C = B*A;
    return
end
if isa(A,'nopvar') || isa(A,'ndopvar')
    A = ndopvar2tensopvar(A);
elseif ~isa(A,'tensopvar') && ~(isa(A,'double') && all(size(A)==1))
    error("Multiplication of tensopvar objects with non-tensopvar objects is currently not supported.")
end
if isa(B,'nopvar') || isa(B,'ndopvar')
    B = ndopvar2tensopvar(B);
elseif ~isa(B,'tensopvar')
    error("Multiplication of tensopvar objects with non-tensopvar objects is currently not supported.")
end

% Check that the inner dimensions of the objects match
Adim = size(A);
Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Check that the spatial domain of the output of B matches the input of A
if isa(A,'tensopvar') 
    if ~isequal(A.depmat2,B.depmat1) || ~isequal(A.dims(:,2),B.dims(:,1))
        error("Function space of input of A must match that of output of B.")
    end
end

% Perform the multiplication
if isa(A,'tensopvar')
    if (size(A.ops,2)==1 && size(B.ops,2)==1) || (A.type(2) && B.type(1))
        % Support 'tensopvar' multiplication only if the operators acts on
        % monomials of degree 1, or if we have proper tensor products, so
        % that we can use the mixed product property
        %   (A o B)(C o D) = AC o BD
        if B.type(1)
            error("Multiplication of tensor-PI operators of alternative type is not supported.")
        end
        C = B;
        C.ops = cell(size(A.ops,1).*size(B.ops,1),size(A.ops,2));
        for k=1:numel(C.ops,1)
            % Take the product of each pair of operators in A and B
            [idx1,idx2,idx3] = ind2sub([size(A.ops,1),size(B.ops,1),size(C.ops,2)],k);
            C.ops{k} = A.ops{idx1,idx3}*B.ops{idx2,idx3};
        end
        C.depmat1 = A.depmat1;
    else
        error("Proposed multiplication of tensor-PI operators is currently not supported.")
    end
else
    % Multiply each of the operators stored in B by A
    C = B;
    for i=1:size(C.ops,1)
        % Multiply the first factor of each term with A
        C.ops{i,1} = A*B.ops{i,1};
    end
end

end