function [logval,idx_out,idx_in] = is_opvar2d(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [logval,idx_out,idx_in] = is_opvar2d(P) checks whether the PI operator P
% has the structure of an 'opvar2d' object,
%   P: R^m0 x L2^mx[x] x L2^my[y] x L2^m2[x,y]
%                           -> R^n0 x L2^nx[x] x L2^ny[y] x L2^n2[x,y].
%
% INPUT
% P:        copvar class object
%
% OUTPUT
% logval:   logical, true if P has at most two spatial variables, and its
%           block rows (columns) are ordered by space: first R, then
%           L2[x], then L2[y], then L2[x,y];
% idx_out:  M x 1 array, the space of each block row: 1 for R, 2 for
%           L2[x], 3 for L2[y], 4 for L2[x,y]. Empty if logval is false;
% idx_in:   N x 1 array, the same for the block columns;
%
% NOTES:
% - x is the variable P.vars{1}, and y is P.vars{2}, the order of the
%   sorted registry. A container with one variable has no y, and its
%   variable is x, whatever its name.
% - Any of the dimensions may be zero.
% - P may have several block rows of each space, e.g. one block row per
%   group of PDE states. Stacked in order, they are the
%   spaces of the 'opvar2d'. A block row of a space before one of an
%   earlier space (e.g. L2[x,y] before L2[x]) would change the order of the
%   components, which an 'opvar2d' cannot represent, so that P is not
%   opvar2d-like. The same holds for the block columns.
% - A block row or column with zero components is ignored, as it holds no
%   component to order.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_opvar2d (copvar)
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/15/2026: Initial coding;

% An 'opvar2d' has two spatial variables, x and y
logval = false;     idx_out = [];       idx_in = [];
if numel(P.vars)>2
    return
end
% Space index 1 + [x in space] + 2*[y in space] of each row (column), with
% the mask padded to two columns for a container with fewer variables.
% Ignoring rows without components, the indices may not decrease.
spaces = {P.space_out,P.space_in};      dims = {P.dim_out,P.dim_in};
idx = cell(1,2);
for k=1:2
    mask = false(size(spaces{k},1),2);
    mask(:,1:size(spaces{k},2)) = spaces{k};
    idx{k} = 1 + mask*[1;2];
    if any(diff(idx{k}(dims{k}(:)>0))<0)
        return
    end
end
[idx_out,idx_in] = idx{:};
logval = true;

end
