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
% MMP, 09/25/2026: The two decision variable lists are reconciled with
%                  'merge_dvar_lists' instead of 'put_on_list': one
%                  unique(...,'stable') gives both operands' row maps, where
%                  'put_on_list' had every block search the union by name
%                  (4.24 of 5.31 s for D + E at q = 4e5 each; the pattern of
%                  every executive, and 'minus' goes through here). Every
%                  block sum is then pointed at the container's one Zd
%                  array, since a block sum carries its own copy of the list.
% MMP, 09/26/2026: Summands over different variable registries are restated
%                  over their union by 'merge_copvar_registry' instead of
%                  raising plus:registryMismatch (an R^n -> R^n operator
%                  such as -gam*Iw has an empty registry, so Hinf_gain with
%                  Tw ~= 0 was refused). Metadata only, flat in q; a
%                  variable on two domains is cdopvar:domConflict.
% MMP, 09/29/2026: A dpvar summand is routed at the top to
%                  'dpvar_op_copvar' with the legacy @dopvar/plus semantics:
%                  a scalar is scalar*I on a square operator (gam*I - Q), a
%                  matrix of the total size a multiplier split by the
%                  dimensions. Before, it reached the dpvar's method or
%                  errored plus:badInput. A numeric summand is still
%                  plus:badInput. minus = plus(A,-B) reaches the branch
%                  unchanged.
% MMP, 09/30/2026: The reconciliation calls 'merge_dvar_pair' (new, private),
%                  the two-operand wrapper of 'merge_dvar_lists' that this
%                  file and 'mtimes' each carried inline. It reads each list
%                  as stored, where A.Zd(:) and B.Zd(:) made up to 5 q-length
%                  copies per sum. Same union, same blocks, same result.
%                  'private/put_on_list' (no caller since 09/26/2026) is
%                  deleted; the '(was)' lines below that call it are history.
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): log entries
%                  condensed, a (was) comment line removed; code unchanged.

% dpvar summand: legacy scalar*I / matrix multiplier, before                % MMP, 09/29/2026
% the container checks below.                                               % MMP, 09/29/2026
if isa(A,'dpvar') || isa(B,'dpvar')                                         % MMP, 09/29/2026
    C = dpvar_op_copvar('plus',A,B);                                        % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026
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
% 'merge_dvar_lists' for why doing it as a by-product of the loop is wrong. % MMP, 09/30/2026
% Merging first also puts every block addition below on its fast lane.
% CA = A.C;   CB = B.C;   Zd = A.Zd(:);                                     % MMP, 09/30/2026 (was)
% if ~isequal(A.Zd(:),B.Zd(:))                                              % MMP, 09/30/2026 (was)
%     Zd = unique([A.Zd(:);B.Zd(:)],'stable');                              % MMP, 09/25/2026 (was)
%     CA = put_on_list(CA,Zd);                                              % MMP, 09/25/2026 (was)
%     CB = put_on_list(CB,Zd);                                              % MMP, 09/25/2026 (was)
% Deleted MMP, 09/30/2026: the rest of this 'if' (09/25/2026: both grids    % MMP, 09/30/2026
% as one, operand 1 = A and 2 = B, one 'merge_dvar_lists' call on           % MMP, 09/30/2026
% {A.Zd(:),B.Zd(:)}) and its 'end'. It is now 'merge_dvar_pair', shared     % MMP, 09/30/2026
% with 'mtimes', which reads no list through a q-length (:) copy.           % MMP, 09/30/2026
[CA,CB,Zd] = merge_dvar_pair(A.C,B.C,A.Zd,B.Zd);                            % MMP, 09/30/2026

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
