function k = cell_of_gamma(G)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% K = CELL_OF_GAMMA(G) returns the linear parameter cell index of each
% multi-index gamma in the rows of G, for a 'sopvar' or 'sdopvar' whose
% parameter cell spans size(G,2) shared variables. The inverse of
% 'gamma_of_cell'.
%
% INPUT
% - G:      r x n3 array, row i a multi-index with entries 1 (multiplier),
%           2 (lower integral) or 3 (upper integral). A 1 x n3 row is one
%           multi-index; r x 0 is r copies of the empty multi-index (n3 = 0);
%
% OUTPUT
% - k:      r x 1 double, k(i) the linear index of cell G(i,:) in the
%           3 x ... x 3 parameter cell array; 1 for the empty multi-index.
%
% NOTES
% Convention (sopvar_implementation_notes.pdf Sec. 1, 3 and 4) as in
% 'gamma_of_cell': MATLAB's column-major linear index, direction 1 fastest,
%
%   k = 1 + sum_t (gamma_t - 1)*3^(t-1),
%
% which is sub2ind([3*ones(1,n3),1],gamma_1,...,gamma_n3) for one row.
% Integer arithmetic, exact for n3 <= 33.
%
% Cost: O(r*n3) on the spatial axis.
%
% See also GAMMA_OF_CELL, SUB2IND.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - cell_of_gamma
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
% Initial coding MMP, 09/30/2026: one public form of the gamma -> cell-index
%                  map, which 'eq_opts_sopvar' (lin_of), 'degbalance_core'
%                  (inline loop) and 'lpi_eq_sdopvar' (sub2ind) each wrote
%                  locally.

% An entry outside 1..3 would address another cell silently; sub2ind
% raises an error there, so this does too.
if any(G(:)<1 | G(:)>3 | G(:)~=fix(G(:)))
    error("Multi-index entries should be 1, 2 or 3.")
end
k = 1 + (G-1)*(3.^(0:size(G,2)-1)).';

end
