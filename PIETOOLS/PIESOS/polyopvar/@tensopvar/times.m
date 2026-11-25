function C = times(A,B)
% C = TIMES(A,B) returns the 'tensopvar' object C representing the
% elementwise product of the multiplier operator defined by A with the
% tensor-PI operator defined by B.
%
% INPUTS
% - A:      m x n 'double' object (or constant 'polynomial') representing
%           some multiplier operator to multiply with B (a 'dpvar' is not
%           supported: the factors cannot depend on decision variables);
% - B:      m x n 'tensopvar' object acting on some degree-d distributed
%           monomial;
%
% OUTPUTS
% - C:      m x n 'tensopvar' object representing the product of the
%           multiplier A with the tensor-PI operator defined by B;
%
% NOTES
% - If A is a 'tensopvar' and B is not, then B.*A is returned. The product
%   of two 'tensopvar' objects is not supported.
% - A multiplies the first factor of each term of B. If A or B has a
%   single row (column), it is repeated to match the number of rows
%   (columns) of the other.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - times
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

% Check that multiplication is supported
if isa(A,'tensopvar') && isa(B,'tensopvar')
    error("Multiplication of tensor-PI operators is not supported.")
elseif isa(A,'tensopvar')
    % Make sure B is the tensor-PI operator
    C = B.*A;
    return
elseif isa(A,'dpvar')
    % The factors of a 'tensopvar' cannot depend on decision variables
    error("Multiplication of a 'tensopvar' with a 'dpvar' is not supported: the factors of a 'tensopvar' cannot depend on decision variables.")
elseif isa(A,'polynomial') && isempty(A.varname)
    A = double(A);
elseif isa(A,'polynomial')
    % The factors ('sopvar') do not support a multiplier in spatial variables
    error("Elementwise multiplication of a 'tensopvar' with a 'polynomial' in spatial variables is not supported.")
elseif ~isa(A,'double')
    error("'tensopvar' objects can only be multiplied with objects of type 'double' or constant 'polynomial'.")
end
if ~B.chkval
    B = initialize(B);
end

% Allow for multiplication of scalar with matrix
Adim = size(A);     Bdim = size(B);
if all(Adim==1)    
    C = B;
    for i=1:size(C.ops,1)
        % Multiply the first factor of each term with A
        C.ops{i,1} = A*B.ops{i,1};
    end
    return
end

% Make sure the row dimensions of the objects match
isdim1 = [false,false];
if Adim(1)~=Bdim(1)
    if Adim(1)==1
        A = ones(Bdim(1),1).*A;
    elseif Bdim(1)==1
        isdim1(1) = true;
    else
        error("Row dimensions of objects to multiply must match.")
    end
end
% Make sure the column dimensions of the objects match
if Adim(2)~=Bdim(2)
    if Adim(2)==1
        A = ones(1,Bdim(2)).*A;
    elseif Bdim(2)==1
        isdim1(2) = true;
    else
        error("Column dimensions of objects to multiply must match.")
    end
end
Adim = size(A);     % A may have been expanded above

% A single row (column) of B is repeated to match the rows (columns) of A
C = B;
if any(isdim1)
    C = repmat(B,1+(Adim(1)-1)*isdim1(1),1+(Adim(2)-1)*isdim1(2));
end
% Multiply the first factor of each term with A, elementwise; a 'sopvar'
% factor as A.*R = sum_j diag(A(:,j))*R*E_jj, E_jj = e_j*e_j'
for i=1:size(C.ops,1)
    R1 = C.ops{i,1};
    if ~isequal(size(R1),Adim)
        error("Elementwise multiplication with a matrix requires the first factor of each term to have the size of the matrix.")
    elseif isa(R1,'sopvar') && isnumeric(A)
        Ejj = zeros(Adim(2));       Ejj(1,1) = 1;
        R1new = diag(A(:,1))*R1*Ejj;
        for j=2:Adim(2)
            Ejj = zeros(Adim(2));   Ejj(j,j) = 1;
            R1new = R1new + diag(A(:,j))*R1*Ejj;
        end
        C.ops{i,1} = R1new;
    else
        C.ops{i,1} = A.*R1;
    end
end

end