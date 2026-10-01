function A = multiindex_grid(vals,order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A = MULTIINDEX_GRID(VALS) returns every tuple obtained by picking one
% entry of VALS{d} in each direction d, one tuple per row, DIRECTION 1
% VARYING FASTEST.
%
% A = MULTIINDEX_GRID(VALS,'first_slowest') returns the same tuples with
% DIRECTION 1 VARYING SLOWEST, the order of kron(Z{1},...,Z{N}).
% 'first_fastest' is the default.
%
% INPUT
% - vals:   1 x N cell of numeric vectors, any orientation;
% - order:  (optional) 'first_fastest' (default) or 'first_slowest';
%
% OUTPUT
% - A:      prod(numel(vals{d})) x N double. N = 0 gives zeros(1,0), the
%           single empty tuple; an empty vals{d} gives zeros(0,N).
%
% NOTES
% The two orders are the two index conventions of the classes:
%
% - 'first_fastest' is the order of the 3 x ... x 3 parameter cell array
%   (sopvar_implementation_notes.pdf Sec. 1 and 3; see 'gamma_of_cell'),
%   and of the basis operators Z_alpha over alpha (Sec. 9.1). With
%   vals{d} = [1,2,3] (or [1,4] in a separable direction) it is the
%   multi-index list of 'possopvar'/'sopquadvar'/'copquadvar', whose linear
%   'include' indices refer to its rows; with vals{d} = 1:3 in every
%   direction, row k is gamma_of_cell(k,N).
% - 'first_slowest' is the order of the monomial basis
%   Z(s) = Z_1(s_1) (x) ... (x) Z_N(s_N) (Sec. 8), first variable outermost.
%   With vals{d} = Z{d} the exponent vectors of a basis, row t is the
%   exponent multi-index of monomial t: the degree table of the basis.
%
% A Kronecker POSITION list over a product of per-variable index sets, in
% first-fastest enumeration order, is multiindex_grid(idx)*st(:) with
% st = kron_strides(nvec).
%
% The entries of A are entries of vals, copied, never computed, so any two
% constructions of the same order agree bit for bit.
%
% Cost: O(prod(numel(vals{d}))*N), on the spatial or monomial axis.
%
% See also GAMMA_OF_CELL, CELL_OF_GAMMA, KRON_STRIDES, KRON_SPLIT, NDGRID.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - multiindex_grid
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
% Initial coding MMP, 09/30/2026: one public product-grid enumeration. It
%                  replaces 'alpha_grid' (copquadvar), the inline copy in
%                  'sopquadvar', 'enum_alpha' and the grid of 'reach'
%                  (eq_opts_sopvar, degbalance_core), and 'degree_table'
%                  (degbalance_core, monomial_gather), written with two
%                  algorithms (repmat, kron) and two orders. The body is
%                  alpha_grid's, run from the other end for 'first_slowest'.
%                  Differences from the copies, only where a copy errors
%                  or is wrong: N = 0 returns zeros(1,0) where enum_alpha
%                  and reach index vals{1} and error; an empty vals{d}
%                  returns zeros(0,N) where alpha_grid errors (0/0
%                  replication) and degbalance_core's degree_table, if a
%                  nonempty basis follows, returns a table with rows and
%                  one column fewer (measured); a ROW exponent vector is
%                  read as a list, where monomial_gather's degree_table
%                  widens the table.

if nargin<2 || strcmp(order,'first_fastest')
    first_fastest = true;
elseif strcmp(order,'first_slowest')
    first_fastest = false;
else
    error("'order' should be 'first_fastest' or 'first_slowest'.")
end

nd = numel(vals);
szv = cellfun(@numel,vals);
nall = prod([szv(:).',1]);
A = zeros(nall,nd);
if nall==0
    return              % an empty value set: no tuple
end
% rep: how many consecutive rows share one value of direction d, i.e. the
% product of the sizes of the directions varying faster than d.
if first_fastest,   dirs = 1:nd;
else,               dirs = nd:-1:1;
end
rep = 1;
for d = dirs
    col = reshape(repmat(reshape(vals{d},1,[]),rep,1),[],1);
    A(:,d) = repmat(col,nall/(rep*szv(d)),1);
    rep = rep*szv(d);
end

end
