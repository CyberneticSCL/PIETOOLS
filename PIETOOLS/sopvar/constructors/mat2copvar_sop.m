function Pop = mat2copvar_sop(Mat,dims,spaces,dom,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = MAT2COPVAR_SOP(MAT,DIMS,SPACES,DOM,OPTIONS) returns the container
% of a spatially constant matrix between given spaces,
%
%   Pop: L_2^{p_1}[s^1] x ... x L_2^{p_N}[s^N]
%                   -> L_2^{q_1}[t^1] x ... x L_2^{q_M}[t^M],
%   (Pop*x)_i = sum_j MAT_ij x_j,
%
% MAT_ij the (i,j) block of MAT, rows cut by q, columns by p. It is the
% container counterpart of 'mat2opvar', with any number of spatial
% variables. A 'double' gives a FIXED multiplier ('copvar'), a 'dpvar' a
% DECISION multiplier ('cdopvar', affine in the dpvar's decision variables),
% as 'mat2opvar' gives opvar or dopvar.
%
% A block is read in the one way the kernel form of Sec. 4 of the sopvar
% document allows a constant: a multiplier (delta) in every variable its two
% spaces share, constant in an output-only variable, and integrated over an
% input-only variable. So R^p -> L2^q[s] is the constant function (Q2),
% L2^p[s] -> L2^q[s] the matrix multiplier (R0), and L2^p[s] -> R^q the
% integral int MAT_ij x(s) ds (Q1), which 'mat2opvar' refuses and
% OPTIONS.mult_only refuses here too.
%
% INPUT
% - Mat:     sum(q) x sum(p) 'double' or 'dpvar', constant in space (a
%            dpvar in a spatial variable is an error);
% - dims:    struct with fields 'out' (M x 1) and 'in' (N x 1), or an M x 1
%            array for the same spaces in and out ('lpivar_cdopvar');
% - spaces:  struct with fields 'out' (1 x M cell) and 'in' (1 x N cell) of
%            cellstr variable names, {} being R^q, or a 1 x M cell for the
%            same spaces in and out ('lpivar_cdopvar', 'poscopvar');
% - dom:     nv x 2 in the order of the SORTED registry of all variables, a
%            1 x 2 row for every variable, or a struct with fields 'vars' and
%            'dom';
% - options: (optional) struct with field
%            - mult_only: true refuses a nonzero integral block (default
%                         false).
%
% OUTPUT
% - Pop:     'copvar' (double) or 'cdopvar' (dpvar), passing 'verify', in
%            canonical multiplier form. Zero blocks are [], except one
%            explicit zero block per row or column that would otherwise
%            have none.
%
% EXAMPLE - Iw, Iz and -gam*Iw of the H-infinity executives:
%   Iw = mat2copvar_sop(eye(nw),nw,{{}},[]);            % R^nw -> R^nw
%   [prog,gam] = lpidecvar(prog,'gam');
%   Gw = mat2copvar_sop(-gam*eye(nw),nw,{{}},[]);      % a cdopvar
%
% NOTES
% Cost: O(M*N) blocks, each O(q_i*p_j) plus the nonzeros of its slice of
% the dpvar's coefficients; O(nnz(Mat.C)) once for a dpvar. Nothing is
% O(number of decision variables) beyond that.
%
% See also EYE_COPVAR_SOP, ZEROS_COPVAR_SOP, MAT2OPVAR, LPIVAR_CDOPVAR,
%          MAT2COPVAR_GRID.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mat2copvar_sop
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
% Initial coding MMP, 09/29/2026. Tier 1c of the container parity map: the
%                constant and decision matrix multipliers over given
%                spaces, replacing opvar2copvar(mat2opvar(...)), which has
%                no N-D route. The work is in 'mat2copvar_grid' (sopvar/
%                misc), which the classes' dpvar branches share.
% MMP, 09/30/2026: Read the spaces with the shared 'parse_copvar_spaces'
%                (sopvar/misc/conventions); 'private/spaces2meta_sop', a
%                verbatim copy of the same parser, is deleted. Same meta on
%                every input; a non-cell space list now raises the parser's
%                message instead of MATLAB's brace-indexing error.

if nargin<5 || isempty(options),    options = struct();     end
mult_only = isfield(options,'mult_only') && options.mult_only;
% meta = spaces2meta_sop(dims,spaces,dom);                                  % MMP, 09/30/2026 (was)
meta = parse_copvar_spaces(dims,spaces,dom);                                % MMP, 09/30/2026
Pop = mat2copvar_grid(Mat,meta,mult_only);

end
