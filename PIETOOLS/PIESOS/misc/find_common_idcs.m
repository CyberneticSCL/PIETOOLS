function [idcs1,idcs2] = find_common_idcs(arr1,arr2,nrange)
% [IDCS1,IDCS2] = FIND_COMMON_IDCS(ARR1,ARR2,NRANGE) determines all indices
% IDCS1 and IDCS2 such that ARR1(IDCS1)=ARR2(IDCS2).
%
% INPUTS
% - arr1, arr2:     m x 1 and n x 1 numeric arrays;
% - nrange:         (optional) scalar; specify it only if arr1 and arr2
%                   take positive integer values in [1,nrange] (e.g.
%                   row/column indices of two conformable
%                   sparse-matrix-triplet structs). The values are then
%                   used directly as labels, instead of being relabeled
%                   by 'unique';
%
% OUTPUTS
% - idcs1,idcs2:    k x 1 numeric arrays for some k<=m*n, such that i is
%                   in idcs1 and j is in idcs2 if and only if
%                       arr1(i) = arr2(j);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - find_common_idcs
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
% DJ, 10/15/2026: Initial coding;

if nargin>=3 && ~isempty(nrange)
    % Caller-supplied range: skip 'unique', build membership directly.
    nnz1 = numel(arr1);         nnz2 = numel(arr2);
    S1 = sparse((1:nnz1)',arr1(:),true,nnz1,nrange);
    S2 = sparse(arr2(:),(1:nnz2)',true,nrange,nnz2);
    [idcs1, idcs2] = find(S1*S2);
    idcs1 = idcs1(:);   idcs2 = idcs2(:);
    return
end

% Assign a common index idx to matching values in arr1 and arr2
[~, ~, idx] = unique([arr1(:); arr2(:)]);
nnz1 = numel(arr1);         nnz2 = numel(arr2);
idx1 = idx(1:nnz1);         idx2 = idx(nnz1+1:end);
m = max(idx);

% Declare SPARSE membership matrices:
%   S1(i,k) = true if arr1(i) has index k
%   S2(k,j) = true if arr2(j) has index k
S1 = sparse((1:nnz1)', idx1, true, nnz1, m);
S2 = sparse(idx2,(1:nnz2)', true, m, nnz2);

% Find all indices for which arr1 and arr2 match
[idcs1, idcs2] = find(S1*S2);
idcs1 = idcs1(:);   idcs2 = idcs2(:);

end