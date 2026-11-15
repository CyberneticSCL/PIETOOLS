function logval = is_opvar1d(P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% logval = is_opvar1d(P) checks whether the PI operator P has the structure
% of a 'dopvar' object,
%   P: R^m x L2^n[a,b] -> R^p x L2^q[a,b],
% so that its fields can be read as those of a 'dopvar' with
% 'getfield_opvar1d'.
%
% INPUT
% P:      cdopvar class object
%
% OUTPUT
% logval: logical, true if P has at most one spatial variable, and every
%         block row (column) of R comes before every block row (column) of
%         L2;
%
% NOTES:
% - m, n, p, q may be zero. A container without a spatial variable maps
%   R^m -> R^p.
% - P may have several block rows of R, and several of L2[a,b], as 'convert'
%   gives one block row per group of PDE states. Stacked in order, they are
%   R^p and L2^q. A block row of R after one of L2 would put components of
%   R^p after those of L2^q, which a 'dopvar' cannot represent, so that P
%   is not opvar-like. The same holds for the block columns.
% - A block row or column with zero components is ignored, as it holds no
%   component to order.
% - Only the metadata of P is read, not its blocks or decision variables.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - is_opvar1d (cdopvar)
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
% DJ, 11/15/2026: Initial coding;

% A 'dopvar' has a single spatial variable, s in [a,b]
logval = false;
if numel(P.vars)>1
    return
end
% With at most one variable, a row (column) is L2[a,b] if its space mask is
% nonempty, and R otherwise. No R row may follow an L2 row; find(...,1) is
% empty if there is no L2 row, and then nothing is checked.
spaces = {P.space_out,P.space_in};      dims = {P.dim_out,P.dim_in};
for k=1:2
    has_comp = dims{k}(:)>0;
    is_L2 = any(spaces{k},2) & has_comp;
    is_R = ~any(spaces{k},2) & has_comp;
    if any(is_R(find(is_L2,1):end))
        return
    end
end
logval = true;

end
