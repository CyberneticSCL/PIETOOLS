function out = monomial_position(map,deg,basis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = MONOMIAL_POSITION(MAP,DEG) returns the 0-based position of each
% exponent in DEG within the one-variable basis that MAP describes, and
% raises an error if the basis does not hold one of them.
%
% OUT = MONOMIAL_POSITION(MAP,DEG,BASIS) names the basis in that error.
%
% INPUT
% - map:    (max exponent+1) x 1 double, the lookup table of the basis:
%           map(p+1) the 0-based position of exponent p, -1 where the basis
%           does not hold p;
% - deg:    array of exponents, any shape;
% - basis:  (optional) string naming the basis in the error message,
%           default "Monomial basis";
%
% OUTPUT
% - out:    numel(deg) x 1 double of 0-based positions.
%
% NOTES
% A missing degree means a basis was enlarged without it, an internal
% inconsistency rather than a user input error, so it raises an error
% rather than returning -1.
%
% Cost: O(numel(deg)), on the coefficient positions being re-indexed.
%
% See also MONOMIAL_POSITION_MAP, KRON_STRIDES.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - monomial_position
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
% Initial coding MMP, 09/30/2026

% Default name of the basis in the error message
if nargin<3
    basis = "Monomial basis";
end
% Look up the exponents in the range of the table; -1 marks any other
deg = double(deg(:));
out = -ones(numel(deg),1);
in_range = deg>=0 & deg==round(deg) & deg+1<=numel(map);
out(in_range) = map(deg(in_range)+1);
% Error on the first exponent the basis does not hold
bad = find(out<0,1);
if ~isempty(bad)
    error(basis+" does not contain degree "+num2str(deg(bad))...
          +"; the bases were not enlarged correctly.")
end

end
