function C = plus(A,B)
% C = PLUS(A,B) returns the 'intop' object C representing the sum of the
% functionals defined by the 'intop' objects A and B.
%
% INPUTS
% - A:      m x n 'intop' object acting on a degree-d distributed
%           monomial;
% - B:      m x n 'intop' object acting on the same degree-d distributed
%           monomial as A;
%
% OUTPUTS
% - C:      m x n 'intop' object representing the sum of the functionals
%           defined by A and B;
%
% NOTES
% - If A or B is the scalar 0 or an empty array, the other is returned.
% - A and B must have the same variables, degrees of the variables
%   (depmat), and domain.
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

% Support addition with 0
if isequal(A,0) || isempty(A)
    C = B;
    return
end
if isequal(B,0) || isempty(B)
    C = A;
    return
end
% Check that the operators are well-defined
if (~A.chkval && ~isvalid(A)) || (~B.chkval && ~isvalid(B))
    error("One of the operators is not properly declared.")
end
% Check that the operators are of the same dimension
if any(size(A)~=size(B))
    error("Matrix dimensions of functionals must match for addition.")
end
% Check that the operators act on a distributed monomial of the same degree
if ~isequal(A.vars,B.vars) || ~isequal(A.depmat,B.depmat)
    error("Functionals must act on distributed monomial of same degree.")
end
% Check that the domain of the spatial variables matches
if ~isequal(A.dom,B.dom)
    error("Domain of integration of the functionals must match.")
end

% Express the functionals in terms of the same integrals
[A,B] = common_intmat(A,B);

% Add the kernels associated with the same integral
C = A;
C.params = cellfun(@(a,b) a+b,A.params,B.params,'UniformOutput',false);
C.chkval = true;

end
