function P = horzcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = HORZCAT(A,B,...) or [A, B, ...] concatenates containers side by side:
% the block columns of A, then those of B, and so on, each block unchanged.
%
%   [A, B] * [x; y] = A*x + B*y,   A: X -> Y,  B: Z -> Y.
%
% INPUTS
% - A, B, ...:  'copvar' objects, or single 'sopvar' blocks, which are taken
%               as 1 x 1 containers. They must agree on every block ROW: the
%               same output space and component count in row i. They may be
%               built over different variable registries; the registries are
%               merged, and a variable given two domains is an error. [] is
%               ignored;
%
% OUTPUTS
% - P:          'copvar' object. If any operand is an 'sdopvar' block, the
%               result is a 'cdopvar' and is built by @cdopvar/horzcat.
%
% NOTES
% Spaces are kept separate, not merged: see 'cat_copvar_grid', which holds
% the logic both container classes share. Without this method MATLAB's
% builtin concatenation silently builds a 1 x K ARRAY of copvar objects,
% which is not an operator and fails only later, elsewhere.
%
% Cost: O(total blocks) plus set operations on the variable registries.
%
% See also VERTCAT, BLKDIAG, CAT_COPVAR_GRID, COPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - horzcat(copvar)
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
% here; promote it and re-dispatch, which reaches @cdopvar/horzcat.
isdec = cellfun(@(a) isa(a,'sdopvar'),varargin);
if any(isdec)
    varargin(isdec) = cellfun(@(a) cdopvar({a}),varargin(isdec),'uni',0);
    P = horzcat(varargin{:});
    return
end

[C,meta] = cat_copvar_grid('h',varargin,'copvar');
P = copvar(C,meta);

end
