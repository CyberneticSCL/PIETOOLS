function C = times(A,B)
% C = TIMES(A,B) returns the 'polyopvar' object C representing the
% elementwise product A.*B of two distributed polynomials, or of a
% constant with a distributed polynomial.
%
% INPUTS
% - A:     m x n 'polyopvar' object, or a constant such that A.*B.C is
%          defined, see @tensopmat/times.m;
% - B:     m x n 'polyopvar' object, or a constant such that A.C.*B is
%          defined, see @tensopmat/times.m;
%
% OUTPUTS
% - C:     'polyopvar' object representing A.*B;
%
% NOTES
% - The monomials of C are the products of each monomial of A with each
%   monomial of B, where the factors of B follow those of A, see
%   'polyopvar'.
% - A.*B is the elementwise product otimes(A,B,1), see @polyopvar/otimes.m,
%   a function of a single copy of the variables. The tensor product, a
%   function of a copy of the variables for each factor, is otimes(A,B,2).
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
% DJ, 10/15/2026: Initial coding;
% DJ, 10/22/2026: Move the body to otimes, calling it with otype=1 for the
%                   elementwise product;

C = otimes(A,B,1);

end
