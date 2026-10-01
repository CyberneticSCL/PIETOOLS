function s = kron_strides(nvec)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S = KRON_STRIDES(NVEC) returns the stride of each variable in the
% monomial index of kron(Z{1},...,Z{N}), numel(Z{d}) = NVEC(d): FIRST
% VARIABLE SLOWEST.
%
% INPUT
% - nvec:   1 x N or N x 1 array, nvec(d) the number of monomials of
%           variable d;
%
% OUTPUT
% - s:      1 x N double, s(d) = prod(nvec(d+1:N)), so that monomial
%           (a_1,...,a_N) (0-based per-variable positions) sits at 0-based
%           position a*s(:) of the composite basis. N = 0 gives ones(1,0).
%
% NOTES
% Convention (sopvar_implementation_notes.pdf Sec. 8): the monomial basis
% of a kernel is the Kronecker product of per-variable bases,
% Z(s) = Z_1(s_1) (x) ... (x) Z_N(s_N), over the variables in the class's
% order (vars.out for ZL, vars.in for ZR), so variable 1 is the outermost
% factor. 'kron_split' inverts the position map.
%
% Cost: O(N).
%
% See also KRON_SPLIT, MULTIINDEX_GRID, KRON.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kron_strides
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
% Initial coding MMP, 09/30/2026: one public form of 'strides_of'
%                  (canonicalize_multiplier, canonical_adjoint_map) and
%                  'strides' (lpivar_cdopvar), identical copies of the
%                  first-variable-slowest rule. Body verbatim.

N = numel(nvec);
s = ones(1,N);
for k = N-1:-1:1
    s(k) = s(k+1)*nvec(k+1);
end

end
