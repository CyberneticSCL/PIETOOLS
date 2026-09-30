function Zop = zeros_copvar_sop(dims,spaces,dom)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZOP = ZEROS_COPVAR_SOP(DIMS,SPACES,DOM) returns the zero operator
%
%   L_2^{p_1}[s^1] x ... x L_2^{p_N}[s^N]
%                   -> L_2^{q_1}[t^1] x ... x L_2^{q_M}[t^M]
%
% as an M x N 'copvar', e.g. the zero D11 or the zero blocks of a
% closed-loop PIE.
%
% INPUT
% - dims, spaces, dom: as for 'mat2copvar_sop' (rectangular allowed).
%
% OUTPUT
% - Zop:    M x N 'copvar' passing 'verify'. Every block is [] except one
%           explicit zero block (degree 0, no content) per row and column
%           that would otherwise have none, placed as 'lpivar_cdopvar'
%           places it: first (i,1) for each row, then (1,j) for each column
%           still empty. A container needs them because 'verify' and
%           copvar(C) read a row's space off a block.
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
% Initial coding MMP, 09/29/2026. Tier 1c: zero operators over given
%                spaces, without opvar2copvar(mat2opvar(...)).

meta = spaces2meta_sop(dims,spaces,dom);
Zop = mat2copvar_grid(sparse(sum(meta.dim_out),sum(meta.dim_in)),meta,false);

end
