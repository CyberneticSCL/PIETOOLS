function Zop = zeros_copvar_sop(dims,spaces,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZOP = ZEROS_COPVAR_SOP(DIMS,SPACES,DOM) returns the zero operator
%
%   L_2^{p_1}[s^1] x ... x L_2^{p_N}[s^N]
%                   -> L_2^{q_1}[t^1] x ... x L_2^{q_M}[t^M]
%
% as an M x N 'copvar'.
%
% INPUT
% - dims:   struct with fields 'out' (M x 1) and 'in' (N x 1), the
%           component counts, or an M x 1 array for the same spaces in and
%           out;
% - spaces: struct with fields 'out' (1 x M cell) and 'in' (1 x N cell) of
%           cellstr variable names, {} being R^q, or a 1 x M cell for the
%           same spaces in and out;
% - dom:    nv x 2 in the order of the SORTED registry of all variables, a
%           1 x 2 row for every variable, or a struct with fields 'vars' and
%           'dom'.
%
% OUTPUT
% - Zop:    M x N 'copvar' with consistent metadata. Every block is []
%           except one explicit zero block (degree 0, no content) per row
%           and column that would otherwise have none: first (i,1) for each
%           row, then (1,j) for each column still empty. A container needs
%           them because the space of a row or column is read off a block.
%
% NOTES
% Cost: O(M*N) block visits and at most M + N - 1 blocks of q_i*p_j
% zeros; no decision variables.
%
% See also MAT2COPVAR_SOP, EYE_COPVAR_SOP.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - zeros_copvar_sop
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
% Initial coding MMP, 09/29/2026

% Container metadata of the spaces, then the zero matrix as its blocks
meta = parse_copvar_spaces(dims,spaces,dom);
Zop = mat2copvar_grid(sparse(sum(meta.dim_out),sum(meta.dim_in)),meta,false);

end
