function C = plus(A,B)
% C = PLUS(A,B) returns the 'polyopvar' object C representing the sum
% of the distributed polynomials by the 'polyopvar' objects A and B
%
% INPUTS
% - A:      m x n 'polyopvar' object
% - B:      m x n 'polyopvar' object
%
% OUTPUS
% - C:      m x n 'polyopvar' object representing the sum of the
%           distributed polynomials defined by A and B
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
% DJ, 01/15/2026: Initial coding

% Make sure the two distributed polynomials are expressed in terms of the
% same basis of distributed monomials.
[A,B] = common_basis(A,B);

% Then, the coefficients defining the sum are given by the sum of the
% coefficients
C = A;
for ii=1:size(C.degmat,1)
    if sum(C.degmat(ii,:))==1
        % If the monomial is just linear, Z_{i}(x) = xk for some k,
        % then we can use the 'nopvar' plus routine to compute the sum of
        % the operators
        %   Ai*xk + Bi*xk = (Ai+Bi)*xk;
        C.C{ii} = A.C{ii} + B.C{ii};
    else
        % If the monomial is more complicated, e.g. Z_{i}(x) = xk*yl,
        % we cannot simply add the coefficient
        %   (Ai{1}*xk)(Ai{2}*yk) + (Bi{1}*xk)(Bi{2}*yk)
        C.C{ii} = [A.C{ii}; B.C{ii}];
    end
end

end