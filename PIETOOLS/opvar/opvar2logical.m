function P_logic = opvar2logical(Pop,tol)
% P_LOGIC = OPVAR2LOGICAL(POP,TOL) returns a 'logical' object specifying
% for each row and column of the input 'opvar' or 'opvar2d' object whether
% it is nonzero.
%
% INPUT
% - Pop:    'opvar' or 'opvar2d' object;
% - tol:    scalar value specifying the tolerance in establishing whether
%           elements of the input operator are nonzero. Defaults to 1e-14;
%
% OUTPUT
% - P_logic:    'logical' array of same size of the input operator,
%               with each row i and column j indicating whether the
%               associated element Pop(i,j) is nonzero up to the specified
%               tolerance.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 07/16/2024
%

% % % Check the inputs.
% Make sure the input object is of appropriate type
if ~isa(Pop,'opvar') && ~isa(Pop,'dopvar') && ~isa(Pop,'opvar2d') && ~isa(Pop,'dopvar2d')
    error("Input operator should be an object of type 'opvar' or 'opvar2d'.")
end
% If no tolerance is specified, set a default value.
if nargin==1
    tol = 1e-14;
end

% % % For each element of the operator, check if it is nonzero.
[nr,nc] = size(Pop);
P_logic = false(nr,nc);
for ii=1:nr
    for jj=1:nc
        P_logic(ii,jj) = ~eq(Pop(ii,jj),0,tol);
    end
end

end