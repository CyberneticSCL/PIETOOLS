function P = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = VERTCAT(A,B,...) or [A; B; ...] stacks containers: the block rows of A,
% then those of B, and so on, each block unchanged.
%
%   [A; B] * x = [A*x; B*x],   A: X -> Y,  B: X -> Z.
%
% INPUTS
% - A, B, ...:  'copvar' objects, or single 'sopvar' blocks, which are taken
%               as 1 x 1 containers. They must agree on every block COLUMN:
%               the same input space and component count in column j. They
%               may be built over different variable registries; the
%               registries are merged, and a variable given two domains is
%               an error. [] is ignored;
%
% OUTPUTS
% - P:          'copvar' object. If any operand is an 'sdopvar' block, the
%               result is a 'cdopvar' and is built by @cdopvar/vertcat.
%
% NOTES
% Spaces are kept separate, not merged: see 'cat_copvar_grid', which holds
% the logic both container classes share. Without this method MATLAB's
% builtin concatenation silently builds a K x 1 ARRAY of copvar objects.
%
% Cost: O(total blocks) plus set operations on the variable registries.
%
% See also HORZCAT, BLKDIAG, CAT_COPVAR_GRID, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - vertcat(copvar)
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

% copvar outranks sdopvar, so a decision BLOCK among the operands lands
% here; promote it and re-dispatch, which reaches @cdopvar/vertcat.
isdec = cellfun(@(a) isa(a,'sdopvar'),varargin);
if any(isdec)
    varargin(isdec) = cellfun(@(a) cdopvar({a}),varargin(isdec),'uni',0);
    P = vertcat(varargin{:});
    return
end

[C,meta] = cat_copvar_grid('v',varargin,'copvar');
P = copvar(C,meta);

end
