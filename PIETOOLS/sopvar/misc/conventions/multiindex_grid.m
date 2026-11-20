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
%   (linear cell index k, direction 1 fastest) and of the basis operators
%   Z_alpha over alpha. With vals{d} = 1:3 in every direction, row k is the
%   multi-index gamma of cell k; with vals{d} = [1,2,3] (or [1,4] in a
%   separable direction), row r is the basis multi-index of block r.
% - 'first_slowest' is the order of the monomial basis
%   Z(s) = Z_1(s_1) (x) ... (x) Z_N(s_N), first variable outermost.
%   With vals{d} = Z{d} the exponent vectors of a basis, row t is the
%   exponent multi-index of monomial t: the degree table of the basis.
%
% A Kronecker POSITION list over a product of per-variable index sets, in
% first-fastest enumeration order, is multiindex_grid(idx)*st(:) with the
% strides st(d) = prod(nvec(d+1:N)).
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
% Initial coding MMP, 09/30/2026

% Order of the directions
if nargin<2 || strcmp(order,'first_fastest')
    first_fastest = true;
elseif strcmp(order,'first_slowest')
    first_fastest = false;
else
    error("'order' should be 'first_fastest' or 'first_slowest'.")
end

% Number of tuples: the product of the sizes of the value sets
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
    % Each value repeated rep times, the pattern tiled over all rows
    col = reshape(repmat(reshape(vals{d},1,[]),rep,1),[],1);
    A(:,d) = repmat(col,nall/(rep*szv(d)),1);
    rep = rep*szv(d);
end

end
