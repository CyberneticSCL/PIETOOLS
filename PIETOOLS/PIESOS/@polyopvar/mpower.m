function C = mpower(A,d)
% C = MPOWER(A,D) returns the 'polyopvar' object C representing the D-th
% power A^D of the scalar-valued distributed polynomial A.
%
% INPUTS
% - A:     1 x 1 'polyopvar' object;
% - d:     nonnegative integer, the order of the power;
%
% OUTPUTS
% - C:     'polyopvar' object representing A^D;
%
% NOTES
% - Only scalar-valued polynomials are supported for now, for which A^D is
%   the elementwise power A.^D, see @polyopvar/power.m. The matrix power of
%   a matrix-valued polynomial is not supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mpower
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
% DJ, 11/06/2026: Initial coding;

if ~isa(A,'polyopvar')
    error("The base of a power of a distributed polynomial must be a 'polyopvar' object.")
elseif ~all(size(A)==1)
    error("Matrix powers of matrix-valued distributed polynomials are currently not supported.")
end
C = A.^d;

end
