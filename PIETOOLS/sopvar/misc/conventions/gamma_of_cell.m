function G = gamma_of_cell(k,n3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% G = GAMMA_OF_CELL(K,N3) returns the multi-index gamma of each linear
% parameter cell index in K, for a 'sopvar' or 'sdopvar' whose parameter
% cell spans N3 shared variables. The inverse is 'cell_of_gamma'.
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
% Convention (sopvar_implementation_notes.pdf Sec. 1, 3 and 4): the
% parameters form a 3 x ... x 3 cell array, one direction per shared
% variable, and cell gamma holds the kernel carrying the indicator I_gamma,
% where the labels 1, 2, 3 are the spec's gamma_k = 0, 1, -1. K is MATLAB's
% column-major linear index into that array, so DIRECTION 1 VARIES FASTEST:
%
%   k = 1 + sum_t (gamma_t - 1)*3^(t-1).
%
% The directions are the shared variables in sorted order, as the classes
% index them (intersect(vars.in,vars.out)).
%
% This is what '[idcs{:}] = ind2sub([3*ones(1,n3),1],k); cell2mat(idcs)'
% computes one cell at a time, with the same values and a 1 x n3 row per
% cell. Here the whole table of a call, gamma_of_cell((1:3^n3)',n3), is one
% vectorized expression. floor((k-1)/3^t) is exact for k < 2^53, so the
% values are those of ind2sub.
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
% Initial coding MMP, 09/30/2026: one public form of the cell-index ->
%                  gamma map, which the canonical-form routines,
%                  'lpivar_cdopvar', 'eq_opts_sopvar' and 'lpi_eq_sdopvar'
%                  each wrote locally (ind2sub + cell2mat, or mod
%                  arithmetic). A copy that drifted would re-label cells
%                  silently in the others.

k = k(:);
% Out-of-range k is refused: the mod arithmetic would wrap it onto another
% cell, and the ind2sub idiom returns a last entry above 3 (measured).
if any(k<1 | k>3^n3 | k~=fix(k))
    error("Parameter cell indices should be integers in 1..3^n3 = "...
          +num2str(3^n3)+".")
end
G = mod(floor((k-1)./3.^(0:n3-1)),3)+1;

end
