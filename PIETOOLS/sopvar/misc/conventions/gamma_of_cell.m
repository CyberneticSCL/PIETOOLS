function G = gamma_of_cell(k,n3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G = GAMMA_OF_CELL(K,N3) returns the multi-index gamma of each linear
% parameter cell index in K, for a 'sopvar' or 'sdopvar' whose parameter
% cell spans N3 shared variables.
%
% INPUT
% - k:      array of linear cell indices, integers in 1..3^n3, any shape;
% - n3:     number of shared variables, the directions of the cell array;
%
% OUTPUT
% - G:      numel(k) x n3 double, row i the multi-index of cell k(i), with
%           entries 1 (multiplier, delta), 2 (lower integral) and 3 (upper
%           integral). A scalar k gives a 1 x n3 row; n3 = 0 gives
%           numel(k) x 0.
%
% NOTES
% Convention: the parameters form a 3 x ... x 3 cell array, one direction
% per shared variable, and cell gamma holds the kernel carrying the
% indicator I_gamma, where the labels 1, 2, 3 stand for gamma_k = 0, 1, -1
% of the indicator. K is the column-major linear index into that array, so
% DIRECTION 1 VARIES FASTEST:
%
%   k = 1 + sum_t (gamma_t - 1)*3^(t-1).
%
% The directions are the shared variables in sorted order,
% intersect(vars.in,vars.out).
%
% The values equal those of
% '[idcs{:}] = ind2sub([3*ones(1,n3),1],k); cell2mat(idcs)', computed for
% all of K at once; floor((k-1)/3^t) is exact for k < 2^53.
%
% Cost: O(numel(k)*n3) on the spatial axis; nothing depends on decision
% variables.
%
% See also CELL_OF_GAMMA, MULTIINDEX_GRID, IND2SUB.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - gamma_of_cell
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

k = k(:);
% Out-of-range k is refused: the mod arithmetic would wrap it onto another
% cell.
if any(k<1 | k>3^n3 | k~=fix(k))
    error("Parameter cell indices should be integers in 1..3^n3 = "...
          +num2str(3^n3)+".")
end
% The base-3 digits of k-1, plus 1, direction 1 the least significant
G = mod(floor((k-1)./3.^(0:n3-1)),3)+1;

end
