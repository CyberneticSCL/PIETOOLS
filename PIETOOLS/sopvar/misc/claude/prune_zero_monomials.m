function [ZLn,ZRn,keepL,keepR] = prune_zero_monomials(ZL,ZR,occL,occR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [ZLN,ZRN,KEEPL,KEEPR] = PRUNE_ZERO_MONOMIALS(ZL,ZR,OCCL,OCCR) drops the
% monomial degrees that carry no coefficient from a pair of tensor bases.
%
% Composition enlarges ZL and ZR to every degree the intermediate products
% can reach, and most of those degrees end up empty: for the 2D heat
% equation at Lyapunov degree 1 the composed bases hold 31 degrees per
% direction while only 13 are occupied, so 32x of the coefficient area is
% padding (sopvar.pdf S4.1 -- the per-direction integrals set the reachable
% degree, not the attained one).
%
% A degree with no stored coefficient represents the zero function, so
% removing it does not change the operator. It does remove one redundant row
% from every equality constraint the operator later enters, which is what
% makes it worth doing: 'sossolve' fails inside SeDuMi's 'pretransfo' on
% highly redundant equality systems.
%
% Pruning is PER DIRECTION because ZL and ZR are tensor bases: the composite
% monomial index is the tensor index over the per-direction lists, so a
% direction-t degree may be dropped only when every composite monomial whose
% t-th index is that degree is empty. Where the occupied set is not a product
% of per-direction sets the prune is partial rather than lossless, which is
% the common case and still recovers most of the padding.
%
% INPUT
% - ZL,ZR: cell arrays of exponent vectors, as stored by 'sopvar';
% - occL:  prod(cellfun(@numel,ZL)) logical, true where some coefficient of
%          some parameter is nonzero at that composite left monomial;
% - occR:  the same over the composite right monomial.
%
% OUTPUT
% - ZLn,ZRn: the bases with the empty degrees removed, in the input shape;
% - keepL:   prod(cellfun(@numel,ZLn)) gather index, new composite monomial
%            j is old composite monomial keepL(j). Lift it over the matrix
%            dimension to reindex a coefficient's rows;
% - keepR:   the same over the columns.
%
% NOTES
% Nothing is pruned, and the gathers are returned as the identity, when
% every degree of every direction is reached -- so a caller on a hot path
% pays one 'find' per parameter and no data movement. Callers should test
% 'numel(keepL)' against the old size rather than reindexing unconditionally.
%
% The canonical multiplier form survives the prune. Content is not moved, so
% a coefficient that was at a permitted (left,right) degree pair is still
% there; in particular right degree 0 in a multiplier direction is occupied
% whenever any multiplier parameter has content, which is exactly when
% 'canonicalize_multiplier' inspects that direction.
%
% See also MONOMIAL_GATHER, CANONICALIZE_MULTIPLIER, SOPVAR, SDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding MMP, 09/11/2026: composition leaves the output bases padded
%                  with degrees that hold no coefficient, at 32x the
%                  coefficient area for a 2D problem, which both inflates
%                  every downstream operation and hands the solver a large
%                  block of redundant equality rows. Shared by
%                  '@sopvar/mtimes' and '@sdopvar/mtimes'.

[ZLn,keepL] = prune_side(ZL,occL);
[ZRn,keepR] = prune_side(ZR,occR);

end


%%
function [Zn,keep] = prune_side(Z,occ)
% Per-direction prune of one tensor basis, with the composite gather.

Zn = Z;
sz = cellfun(@numel,Z);     sz = sz(:).';
NZ = prod([sz,1]);
keep = (1:NZ).';

% Nothing stored at all, or nothing empty: leave the basis alone. Refusing
% the all-empty case keeps a zero operator on its declared basis rather than
% collapsing it to a shape no caller expects.
if isempty(sz) || ~any(occ) || all(occ)
    return
end

% Decode the composite index. The monomial vector is kron(Z{1},...,Z{N})
% with the FIRST direction outermost, so direction t has stride
% prod(sz(t+1:end)) -- the same convention as 'monomial_gather'.
N = numel(sz);
rem = find(occ(:)) - 1;
kd  = cell(1,N);
for t = 1:N
    str   = prod(sz(t+1:end));
    kd{t} = unique(floor(rem/str)) + 1;
    rem   = mod(rem,str);
end

if isequal(cellfun(@numel,kd),sz)
    return                      % every degree of every direction is reached
end

% Rebuild the composite gather over the kept degrees, first direction
% outermost so the new index has the same tensor layout as the old one.
old = 0;
for t = 1:N
    str = prod(sz(t+1:end));
    old = kron(old,ones(numel(kd{t}),1)) + ...
          kron(ones(numel(old),1),(kd{t}(:)-1)*str);
    Zn{t} = Z{t}(kd{t});
    Zn{t} = Zn{t}(:);           % a basis must be a column; see 'plus'
end
keep = old + 1;

end
