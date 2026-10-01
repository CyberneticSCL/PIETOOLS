function map = monomial_position_map(Zvec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAP = MONOMIAL_POSITION_MAP(ZVEC) returns a lookup table from degree to
% 0-based position in a one-variable monomial basis, for
% 'monomial_position'.
%
% INPUT
% - Zvec:   vector of distinct nonnegative integer exponents, the basis
%           Z_d(s_d) = s_d.^Zvec of one variable;
%
% OUTPUT
% - map:    (max(Zvec)+1) x 1 double, map(p+1) the 0-based position of
%           exponent p in Zvec, and -1 where Zvec does not hold p. An empty
%           Zvec gives zeros(0,1).
%
% NOTES
% A direct table over the shifted degree: the exponents of one variable
% are few and small, so it costs O(max(Zvec)) and each lookup O(1).
%
% Cost: O(max(Zvec)+numel(Zvec)), on the monomial axis.
%
% See also MONOMIAL_POSITION, KRON_STRIDES.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - monomial_position_map
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
% Initial coding MMP, 09/30/2026: one public form of 'degree_lookup', whose
%                  two copies had diverged: 'canonical_adjoint_map' rejects
%                  a negative or non-integer exponent, 'canonicalize_
%                  multiplier' does not and then fails inside MATLAB
%                  indexing on the same input. Body is the rejecting copy,
%                  verbatim; the outputs of the two agree on every basis
%                  either accepts.

d = double(Zvec(:));
if isempty(d)
    map = zeros(0,1);
    return
end
if any(d<0) || any(d~=round(d))
    error("Monomial exponents should be nonnegative integers.")
end
map = -ones(max(d)+1,1);
map(d+1) = 0:numel(d)-1;

end
