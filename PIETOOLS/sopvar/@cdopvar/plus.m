function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = PLUS(A,B) adds two containers blockwise, (A+B)_{ij} = A_{ij} + B_{ij}.
%
% INPUTS
% - A, B:   'cdopvar' objects, or one of them a 'copvar', on the same grid,
%           spaces and dimensions. A 'copvar' operand is promoted, so a
%           decision operator can be added to a fixed one. A numeric or
%           'dpvar' scalar c is c*I on a square operator (gam*I - Q); a
%           'dpvar' matrix of the total size is a multiplier split by the
%           dimensions;
% OUTPUTS
% - C:      'cdopvar' object equal to A + B;
%
% NOTES
% Zero blocks pass through: [] + X is X and [] + [] stays []. A zero block
% carries no basis, so there is nothing to align. This is why the container
% owns the space metadata: a sum of two absent blocks would otherwise have
% no way to say what it maps between. It is also how a MIXED container
% arises: adding a fixed operator to a decision operator with a zero block
% leaves a 'sopvar' block among 'sdopvar' ones, which is legal here.
%
% Spaces are compared as masks over the shared registry, so the canonical
% per-block order of vars.out and vars.in does not enter.
%
% Summands over different registries are first restated over their sorted
% union; an R^n -> R^n summand such as -gam*Iw has an empty registry.
% Metadata only, O(nv).
%
% Cost: M*N block additions, plus one decision variable reconciliation for
% the whole sum rather than one per block. Each block addition unions the
% two bases for that block only, so nothing is inflated to a row or column
% union. Block sums stay pairwise: a batch sum gains only for three or more
% operands sharing one synchronization.
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

% dpvar summand: scalar*I or a matrix multiplier, before
% the container checks below.
% A numeric scalar is the same as a constant dpvar scalar: scalar*I.
if isa(A,'dpvar') || isa(B,'dpvar') || ...
        (isnumeric(A) && isscalar(A)) || (isnumeric(B) && isscalar(B))
    C = dpvar_op_copvar('plus',A,B);
    return
end
% A 'copvar' operand is promoted, so the rest of this routine sees one type.
% Promotion copies the block grid and metadata and leaves Zd empty.
if isa(A,'copvar'),     A = cdopvar(A);     end
if isa(B,'copvar'),     B = cdopvar(B);     end
if ~isa(A,'cdopvar') || ~isa(B,'cdopvar')
    error('plus:badInput','Summands must be cdopvar or copvar objects.')
end
if ~isequal(size(A,'block'),size(B,'block'))
    error('plus:gridMismatch','Summands have block grids %s and %s.',...
        mat2str(size(A,'block')),mat2str(size(B,'block')))
end
% One registry for both, so the masks below compare by name; a variable
% on two domains is an error.
[A,B] = merge_copvar_registry('cdopvar','plus',A,B);
if ~isequal(A.space_out,B.space_out) || ~isequal(A.space_in,B.space_in)
    error('plus:spaceMismatch','Summands map between different spaces.')
end
if ~isequal(A.dim_out(:),B.dim_out(:)) || ~isequal(A.dim_in(:),B.dim_in(:))
    error('plus:dimMismatch','Summands have different component dimensions.')
end

% Reconcile the decision variables once, before the block loop: a block
% whose partner is absent would otherwise keep the list of its operand.
% Merging first also puts every block addition below on its fast lane.
[CA,CB,Zd] = merge_dvar_pair(A.C,B.C,A.Zd,B.Zd);

Cc = cell(size(CA));
for ii = 1:numel(Cc)
    if isempty(CA{ii})
        Cc{ii} = CB{ii};
    elseif isempty(CB{ii})
        Cc{ii} = CA{ii};
    else
        Cc{ii} = CA{ii} + CB{ii};
    end
    % A block sum carries its own copy of the list; share the one array.
    if isa(Cc{ii},'sdopvar'),   Cc{ii}.Zd = Zd;     end
end

% The sum is on the same spaces and dimensions as its summands, so the
% metadata is passed through rather than re-derived.
meta = metadata(A);     meta.Zd = Zd;
C = cdopvar(Cc,meta);

end
