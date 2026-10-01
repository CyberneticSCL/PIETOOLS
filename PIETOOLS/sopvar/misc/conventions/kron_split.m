function a = kron_split(idx,stride)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A = KRON_SPLIT(IDX,STRIDE) splits 0-based positions in the monomial basis
% kron(Z{1},...,Z{N}) into 0-based per-variable positions.
%
% INPUT
% - idx:    array of 0-based composite positions, integers in
%           0..prod(nvec)-1, any shape;
% - stride: 1 x N strides, kron_strides(nvec);
%
% OUTPUT
% - a:      numel(idx) x N double, a(i,d) the 0-based position in Z{d} of
%           composite monomial idx(i), so idx(i) = a(i,:)*stride(:).
%
% NOTES
% Convention (sopvar_implementation_notes.pdf Sec. 8): first variable
% slowest; see 'kron_strides'. Integer arithmetic on doubles, exact below
% 2^53. A position outside 0..prod(nvec)-1 is not detected: column 1 then
% exceeds nvec(1)-1.
%
% Cost: O(numel(idx)*N).
%
% See also KRON_STRIDES, MULTIINDEX_GRID.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kron_split
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
% Initial coding MMP, 09/30/2026: one public form of 'split_index'
%                  (canonicalize_multiplier, canonical_adjoint_map),
%                  identical copies. Body verbatim.

N = numel(stride);
a = zeros(numel(idx),N);
rem = idx(:);
for k = 1:N
    a(:,k) = floor(rem/stride(k));
    rem = rem - a(:,k)*stride(k);
end

end
