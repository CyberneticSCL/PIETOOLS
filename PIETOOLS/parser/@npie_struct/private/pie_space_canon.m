function C = pie_space_canon(depmat,dims)
% C = PIE_SPACE_CANON(DEPMAT,DIMS) returns a canonical description of a
% direct sum of function spaces, so that different orderings and groupings
% of the same spaces can be compared using isequal.
%
% INPUT
% - depmat: M x N logical array, with element (i,k) true if space i depends
%           on spatial variable k;
% - dims:   M x 1 array of integers, specifying the number of components in
%           each space;
%
% OUTPUT
% - C:      K x (N+1) array, each row of which provides the depmat and the
%           total number of components of all spaces with this depmat. The
%           rows are sorted by depmat, and spaces without any components
%           are omitted.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pie_space_canon
%
% Copyright (C)2026  PIETOOLS Team
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
% DJ, 10/29/2026: Initial coding;

dims = dims(:);
nvars = size(depmat,2);
keep = dims>0;
depmat = double(depmat(keep,:));
dims = dims(keep);
if isempty(dims)
    C = zeros(0,nvars+1);
    return
elseif nvars==0
    C = sum(dims);
    return
end
[depmat_uni,~,idx] = unique(depmat,'rows');
C = [depmat_uni,accumarray(idx,dims)];

end