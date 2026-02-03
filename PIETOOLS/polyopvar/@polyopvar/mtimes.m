function C = mtimes(A,B)
% C = mtimes(A,B) returns the 'polyopvar' object C representing the product
% of two polyopvar objects.
%
% INPUTS
% - A:     polyopvar object; 
% - B:     polyopvar object;
%
% OUTPUTS
% - C:     polyopvar object representing A*B.
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
% CR, 01/28/2026: Initial coding
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% express A, B in terms of a common monomial basis.
[A,B] = common_basis(A,B);
C = polyopvar();
C.varname = A.varname;
C.pvarname = A.pvarname;
C.dom = A.dom; % assumes function spaces of A and B are the same.
C.varmat = A.varmat;
C.degmat = repelem(A.degmat, size(B.degmat,1), 1) + repmat(B.degmat,size(A.degmat,1), 1); % REMOVE DUPLICATE MONOMIALS? Assumed to be present for tensopvar_mtimes.
C.C = tensopvar_mtimes(A,B);

end