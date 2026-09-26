function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = PLUS(A,B) adds two 'copvar' objects blockwise, (A+B)_{ij} = A_{ij} +
% B_{ij}, following Sec. 8.2 of the sopvar document.
%
% INPUTS
% - A, B:   'copvar' objects on the same grid, spaces and dimensions;
% OUTPUTS
% - C:      'copvar' object equal to A + B;
%
% NOTES
% Zero blocks pass through: [] + X is X and [] + [] stays []. A zero block
% carries no basis, so there is nothing to align. This is why the container
% owns the space metadata - a sum of two absent blocks would otherwise have
% no way to say what it maps between.
%
% Spaces are compared as masks over the shared registry, so the canonical
% per-block order of vars.out and vars.in does not enter.
%
% Summands over different registries are first restated over their sorted   % MMP, 09/26/2026
% union, as concatenation does ('merge_copvar_registry'): an R^n -> R^n     % MMP, 09/26/2026
% summand such as -gam*Iw has an empty registry. Metadata only, O(nv).      % MMP, 09/26/2026
%
% Cost: M*N block additions. Each unions the two bases for that block only,
% so nothing is inflated to a row or column union; see the Sec. 8.3.1 note
% in 'copvar'. Measured over the subset lattice at 2 components and degree
% 1, the whole loop costs 0.032 s at 4 spaces and 0.062 s at 8, so the
% Sec. 8.2 for-loop is not a bottleneck here - unlike composition. Block
% sums stay pairwise. A decision container reconciles decision variables
% here as well; see @cdopvar/plus.
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
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Moved from
%                  @mopvar/ with the class.
% MMP, 09/26/2026: Summands over different variable registries are restated
%                  over their union by 'merge_copvar_registry' instead of
%                  raising plus:registryMismatch. An R^n -> R^n operator
%                  (-gam*Iw) has an empty registry and a product such as
%                  PB'*Tw keeps its left factor's, so compatible summands
%                  were refused, blocking Hinf_gain with Tw ~= 0. The
%                  registry is metadata: blocks untouched, cost O(nv), flat
%                  in q. A variable on two domains is copvar:domConflict.

if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('plus:badInput','Both summands must be copvar objects.')
end
if ~isequal(size(A),size(B))
    error('plus:gridMismatch','Summands have block grids %s and %s.',...
        mat2str(size(A)),mat2str(size(B)))
end
% if ~isequal(A.vars,B.vars) || ~isequal(A.dom,B.dom)                       % MMP, 09/26/2026 (was)
%     error('plus:registryMismatch',['Summands are built on different variable '...
%         'registries; rebuild them over a common set of variables and domains.']) % MMP, 09/26/2026 (was)
% end                                                                       % MMP, 09/26/2026 (was)
% One registry for both, so the masks below compare by name; the domain     % MMP, 09/26/2026
% conflict check lives there too.                                           % MMP, 09/26/2026
[A,B] = merge_copvar_registry('copvar','plus',A,B);                         % MMP, 09/26/2026
if ~isequal(A.space_out,B.space_out) || ~isequal(A.space_in,B.space_in)
    error('plus:spaceMismatch','Summands map between different spaces.')
end
if ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:))
    error('plus:dimMismatch','Summands have different component dimensions.')
end

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
