function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = PLUS(A,B) adds two containers blockwise, (A+B)_{ij} = A_{ij} + B_{ij},
% following Sec. 8.2 of the sopvar document.
%
% INPUTS
% - A, B:   'cdopvar' objects, or one of them a 'copvar', on the same grid,
%           spaces and dimensions. A 'copvar' operand is promoted, so a
%           decision operator can be added to a fixed one;
% OUTPUTS
% - C:      'cdopvar' object equal to A + B;
%
% NOTES
% Zero blocks pass through: [] + X is X and [] + [] stays []. A zero block
% carries no basis, so there is nothing to align. This is why the container
% owns the space metadata - a sum of two absent blocks would otherwise have
% no way to say what it maps between. It is also how a MIXED container
% arises: adding a fixed operator to a decision operator with a zero block
% leaves a 'sopvar' block among 'sdopvar' ones, which is legal here.
%
% Spaces are compared as masks over the shared registry, so the canonical
% per-block order of vars.out and vars.in does not enter.
%
% Summands over different registries are first restated over their sorted   % MMP, 09/26/2026
% union, as concatenation does ('merge_copvar_registry'): an R^n -> R^n     % MMP, 09/26/2026
% summand such as -gam*Iw has an empty registry. Metadata only, O(nv).      % MMP, 09/26/2026
%
% Cost: M*N block additions, plus one decision variable reconciliation for
% the whole sum rather than one per block. Each block addition unions the
% two bases for that block only, so nothing is inflated to a row or column
% union; see the Sec. 8.3.1 note in 'copvar'. Block sums stay pairwise:
% 'plus_batch' wins only for three or more operands sharing one
% synchronization, which is the inner sum of 'mtimes'.
%
% See also MTIMES, CTRANSPOSE, CDOPVAR, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - plus(cdopvar)
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
%                  @mdopvar/ with the class.
% MMP, 09/25/2026: Reconcile the two decision variable lists with
%                  'merge_dvar_lists' instead of 'put_on_list'. One
%                  unique(...,'stable') then gives both operands' row maps,
%                  where 'put_on_list' had every block search the union by
%                  name: measured on two separately declared containers at
%                  q = 4e5 each, D + E took 5.31 s, 4.24 s of it that per-block
%                  'ismember'. This is the pattern of every executive
%                  (Top'*Qop - Pop, Dop + Deop), and 'minus' goes through
%                  here. Also point every block sum at the container's one Zd
%                  array afterwards: a block sum carries a fresh copy of the
%                  list, one per block, where the container invariant is one.
% MMP, 09/26/2026: Summands over different variable registries are restated
%                  over their union by 'merge_copvar_registry' instead of
%                  raising plus:registryMismatch. An R^n -> R^n operator
%                  (-gam*Iw) has an empty registry and a product such as
%                  PB'*Tw keeps its left factor's, so compatible summands
%                  were refused, blocking Hinf_gain with Tw ~= 0. The
%                  registry is metadata: blocks and Zd untouched, cost
%                  O(nv), flat in q. A variable on two domains is
%                  cdopvar:domConflict.

% A 'copvar' operand is promoted, so the rest of this routine sees one type.
% Promotion copies the block grid and metadata and leaves Zd empty.
if isa(A,'copvar'),     A = cdopvar(A);     end
if isa(B,'copvar'),     B = cdopvar(B);     end
if ~isa(A,'cdopvar') || ~isa(B,'cdopvar')
    error('plus:badInput','Summands must be cdopvar or copvar objects.')
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
[A,B] = merge_copvar_registry('cdopvar','plus',A,B);                        % MMP, 09/26/2026
if ~isequal(A.space_out,B.space_out) || ~isequal(A.space_in,B.space_in)
    error('plus:spaceMismatch','Summands map between different spaces.')
end
if ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:))
    error('plus:dimMismatch','Summands have different component dimensions.')
end

% Reconcile the decision variables once, before the block loop; see
% 'put_on_list' for why doing it as a by-product of the loop is wrong.
% Merging first also puts every block addition below on its fast lane.
CA = A.C;   CB = B.C;   Zd = A.Zd(:);
if ~isequal(A.Zd(:),B.Zd(:))
%     Zd = unique([A.Zd(:);B.Zd(:)],'stable');                              % MMP, 09/25/2026 (was)
%     CA = put_on_list(CA,Zd);                                              % MMP, 09/25/2026 (was)
%     CB = put_on_list(CB,Zd);                                              % MMP, 09/25/2026 (was)
    % Both grids as one, operand 1 = A and 2 = B, so the one sort in        % MMP, 09/25/2026
    % 'merge_dvar_lists' yields both row maps.                              % MMP, 09/25/2026
    nA = numel(CA);                                                         % MMP, 09/25/2026
    src = [ones(1,nA), 2*ones(1,numel(CB))];                                % MMP, 09/25/2026
    [Cab,Zd] = merge_dvar_lists([CA(:);CB(:)]',{A.Zd(:),B.Zd(:)},src);      % MMP, 09/25/2026
    CA = reshape(Cab(1:nA),size(CA));                                       % MMP, 09/25/2026
    CB = reshape(Cab(nA+1:end),size(CB));                                   % MMP, 09/25/2026
end

Cc = cell(size(CA));
for ii = 1:numel(Cc)
    if isempty(CA{ii})
        Cc{ii} = CB{ii};
    elseif isempty(CB{ii})
        Cc{ii} = CA{ii};
    else
        Cc{ii} = CA{ii} + CB{ii};
    end
    % A block sum carries its own copy of the list; share the one array.    % MMP, 09/25/2026
    if isa(Cc{ii},'sdopvar'),   Cc{ii}.Zd = Zd;     end                     % MMP, 09/25/2026
end

% The sum is on the same spaces and dimensions as its summands, so the
% metadata is passed through rather than re-derived.
meta = metadata(A);     meta.Zd = Zd;
C = cdopvar(Cc,meta);

end
