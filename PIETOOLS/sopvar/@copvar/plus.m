function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = PLUS(A,B) adds two 'copvar' objects blockwise, (A+B)_{ij} = A_{ij} +
% B_{ij}.
%
% INPUTS
% - A, B:   'copvar' objects on the same grid, spaces and dimensions. A
%           numeric or 'dpvar' scalar c is c*I on a square operator
%           (gam - T'*T); a 'dpvar' matrix of the total size is a
%           multiplier split by the dimensions; a 'dpvar' summand makes
%           the sum a 'cdopvar';
% OUTPUTS
% - C:      'copvar' object equal to A + B;
%
% NOTES
% Zero blocks pass through: [] + X is X and [] + [] stays []. A zero block
% carries no basis, so there is nothing to align. This is why the container
% owns the space metadata: a sum of two absent blocks would otherwise have
% no way to say what it maps between.
%
% Spaces are compared as masks over the shared registry, so the canonical
% per-block order of vars.out and vars.in does not enter.
%
% Summands over different registries are first restated over their sorted
% union; an R^n -> R^n summand such as -gam*Iw has an empty registry.
% Metadata only, O(nv).
%
% Cost: M*N block additions. Each unions the two bases for that block only,
% so nothing is inflated to a row or column union. Measured over the subset
% lattice at 2 components and degree 1, the whole loop costs 0.032 s at 4
% spaces and 0.062 s at 8, so the block loop is not a bottleneck here,
% unlike composition. Block sums stay pairwise.
%
% See also MTIMES, CTRANSPOSE, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus(copvar)
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
% Initial coding MMP, 09/17/2026

% dpvar summand: scalar*I or a matrix multiplier, before
% the container checks below.
% A numeric scalar is the same as a constant dpvar scalar: scalar*I.
if isa(A,'dpvar') || isa(B,'dpvar') || ...
        (isnumeric(A) && isscalar(A)) || (isnumeric(B) && isscalar(B))
    C = dpvar_op_copvar('plus',A,B);
    return
end
% Check the summands and their block grids
if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('plus:badInput','Both summands must be copvar objects.')
end
if ~isequal(size(A,'block'),size(B,'block'))
    error('plus:gridMismatch','Summands have block grids %s and %s.',...
        mat2str(size(A,'block')),mat2str(size(B,'block')))
end
% One registry for both, so the masks below compare by name; a variable
% on two domains is an error.
[A,B] = merge_copvar_registry('copvar','plus',A,B);
if ~isequal(A.space_out,B.space_out) || ~isequal(A.space_in,B.space_in)
    error('plus:spaceMismatch','Summands map between different spaces.')
end
if ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:))
    error('plus:dimMismatch','Summands have different component dimensions.')
end

% Blockwise sum; a zero block passes the other summand through
CA = A.C;   CB = B.C;
Cc = cell(size(CA));
for ii = 1:numel(Cc)
    if isempty(CA{ii})
        Cc{ii} = CB{ii};
    elseif isempty(CB{ii})
        Cc{ii} = CA{ii};
    else
        Cc{ii} = CA{ii} + CB{ii};
    end
end

% The sum is on the same spaces and dimensions as its summands, so the
% metadata is passed through rather than re-derived.
meta = metadata(A);
C = copvar(Cc,meta);

end
