function P = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = VERTCAT(A,B,...) or [A; B; ...] stacks containers: the block rows of A,
% then those of B, and so on, each block unchanged.
%
%   [A; B] * x = [A*x; B*x],   A: X -> Y,  B: X -> Z.
%
% INPUTS
% - A, B, ...:  'cdopvar' or 'copvar' objects, or single 'sopvar'/'sdopvar'
%               blocks, which are taken as 1 x 1 containers. They must agree
%               on every block COLUMN: the same input space and component
%               count in column j. They may be built over different variable
%               registries; the registries are merged, and a variable given
%               two domains is an error. [] is ignored;
%
% OUTPUTS
% - P:          'cdopvar' object, every block on one decision variable list.
%
% NOTES
% The grid and metadata come from 'cat_copvar_grid', shared with 'copvar';
% the decision variable lists are reconciled by 'merge_dvar_lists', the rule
% '@cdopvar/plus' applies.
%
% Cost: O(total blocks) when the operands share a decision variable list,
% which is the normal case. Otherwise one stable union, O(q log q), gives
% every operand's row map, and each block moves in O(nnz(B)) with no
% per-block name search; see 'merge_dvar_lists'. Measured 0.99 s at
% q = 4e5 per operand.
%
% See also HORZCAT, BLKDIAG, CAT_COPVAR_GRID, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - vertcat(cdopvar)
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
% Initial coding MMP, 09/25/2026

[C,meta,Zds,src] = cat_copvar_grid('v',varargin,'cdopvar');
[C,meta.Zd] = merge_dvar_lists(C,Zds,src);
P = cdopvar(C,meta);

end
