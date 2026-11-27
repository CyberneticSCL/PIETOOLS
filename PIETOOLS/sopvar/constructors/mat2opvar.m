function Pop = mat2opvar(Mat,dims,spaces,dom,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POP = MAT2OPVAR(MAT,DIMS,SPACES,DOM,OPTIONS) returns the container
% of a spatially constant matrix between given spaces,
%
%   Pop: L_2^{p_1}[s^1] x ... x L_2^{p_N}[s^N]
%                   -> L_2^{q_1}[t^1] x ... x L_2^{q_M}[t^M],
%   (Pop*x)_i = sum_j MAT_ij x_j,
%
% MAT_ij the (i,j) block of MAT, rows cut by q, columns by p, with any
% number of spatial variables. A 'double' gives a FIXED multiplier
% ('opvar'), a 'dpvar' a DECISION multiplier ('dopvar', affine in the
% decision variables of the dpvar).
%
% A block is read in the one way the kernel form allows a constant: a
% multiplier (delta) in every variable its two spaces share, constant in an
% output-only variable, and integrated over an input-only variable. So
% R^p -> L2^q[s] is the constant function (Q2), L2^p[s] -> L2^q[s] the
% matrix multiplier (R0), and L2^p[s] -> R^q the integral
% int MAT_ij x(s) ds (Q1), which OPTIONS.mult_only refuses.
%
% INPUT
% - Mat:     sum(q) x sum(p) 'double' or 'dpvar', constant in space (a
%            dpvar in a spatial variable is an error);
% - dims:    struct with fields 'out' (M x 1) and 'in' (N x 1), or an M x 1
%            array for the same spaces in and out;
% - spaces:  struct with fields 'out' (1 x M cell) and 'in' (1 x N cell) of
%            cellstr variable names, {} being R^q, or a 1 x M cell for the
%            same spaces in and out;
% - dom:     nv x 2 in the order of the SORTED registry of all variables, a
%            1 x 2 row for every variable, or a struct with fields 'vars' and
%            'dom';
% - options: (optional) struct with field
%            - mult_only: true refuses a nonzero integral block (default
%                         false).
%
% OUTPUT
% - Pop:     'opvar' (double) or 'dopvar' (dpvar), with consistent
%            metadata, in canonical multiplier form. Zero blocks are [],
%            except one explicit zero block per row or column that would
%            otherwise have none.
%
% EXAMPLE - the identity and a decision multiplier on R^nw:
%   Iw = mat2opvar(eye(nw),nw,{{}},[]);            % R^nw -> R^nw
%   [prog,gam] = lpidecvar(prog,'gam');
%   Gw = mat2opvar(-gam*eye(nw),nw,{{}},[]);      % a dopvar
%
% NOTES
% Cost: O(M*N) blocks, each O(q_i*p_j) plus the nonzeros of its slice of
% the coefficients of the dpvar; O(nnz(Mat.C)) once for a dpvar. Nothing is
% O(number of decision variables) beyond that.
%
% See also EYE_OPVAR, ZEROS_OPVAR, MAT2OPVAR, LPIVAR,
%          MAT2OPVAR_GRID.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mat2opvar
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

% Options; the container metadata of the spaces; the blocks of Mat
if nargin<5 || isempty(options),    options = struct();     end
mult_only = isfield(options,'mult_only') && options.mult_only;
meta = parse_opvar_spaces(dims,spaces,dom);
Pop = mat2opvar_grid(Mat,meta,mult_only);

end
