function out = monomial_position(map,deg,basis)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUT = MONOMIAL_POSITION(MAP,DEG) returns the 0-based position of each
% exponent in DEG within the one-variable basis that MAP describes, and
% raises an error if the basis does not hold one of them.
%
% OUT = MONOMIAL_POSITION(MAP,DEG,BASIS) names the basis in that error.
%
% INPUT
% - map:    lookup table from 'monomial_position_map';
% - deg:    array of exponents, any shape;
% - basis:  (optional) string naming the basis in the error message,
%           default "Monomial basis";
%
% OUTPUT
% - out:    numel(deg) x 1 double of 0-based positions.
%
% NOTES
% The error means a caller's basis enlargement missed a degree, an internal
% inconsistency rather than a user input error, so it is raised rather than
% returned as -1.
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
% Initial coding MMP, 09/30/2026: one public form of 'lookup'
%                  (canonicalize_multiplier, canonical_adjoint_map), whose
%                  copies differ only in the first words of the error
%                  message; BASIS carries that difference. Not named
%                  'lookup': which -all finds three toolbox class methods
%                  of that name, and it does not say what is looked up.
%                  Body verbatim.

if nargin<3
    basis = "Monomial basis";
end
deg = double(deg(:));
out = -ones(numel(deg),1);
in_range = deg>=0 & deg==round(deg) & deg+1<=numel(map);
out(in_range) = map(deg(in_range)+1);
bad = find(out<0,1);
if ~isempty(bad)
    error(basis+" does not contain degree "+num2str(deg(bad))...
          +"; the bases were not enlarged correctly.")
end

end
