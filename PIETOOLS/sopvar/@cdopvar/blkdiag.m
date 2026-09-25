function P = blkdiag(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = BLKDIAG(A,B,...) places containers on the block diagonal, with
% structurally zero blocks off it:
%
%   blkdiag(A,B) * [x; y] = [A*x; B*y].
%
% INPUTS
% - A, B, ...:  'cdopvar' or 'copvar' objects, or single 'sopvar'/'sdopvar'
%               blocks, which are taken as 1 x 1 containers. No agreement
%               between them is needed. Their variable registries are
%               merged, and a variable given two domains is an error. [] is
%               ignored;
%
% OUTPUTS
% - P:          'cdopvar' object, every block on one decision variable list.
%
% NOTES
% The off-diagonal blocks are [], the container's zero block. The grid and
% metadata come from 'cat_copvar_grid', shared with 'copvar'; the decision
% variable lists are reconciled by 'merge_dvar_lists', the rule
% '@cdopvar/plus' applies. Operands declared by separate calls - two
% Lyapunov operators, say - have different lists, so here the union is the
% common case rather than the exception.
%
% Cost: O(total blocks) when the operands share a decision variable list.
% Otherwise one stable union, O(q log q), gives every operand's row map,
% and each block moves in O(nnz(B)); see 'merge_dvar_lists'.
%
% See also HORZCAT, VERTCAT, CAT_COPVAR_GRID, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - blkdiag(cdopvar)
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

[C,meta,Zds,src] = cat_copvar_grid('d',varargin,'cdopvar');
[C,meta.Zd] = merge_dvar_lists(C,Zds,src);
P = cdopvar(C,meta);

end
