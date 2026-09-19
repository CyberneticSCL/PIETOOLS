function At = ctranspose(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At = CTRANSPOSE(A) returns the adjoint of a 'mdopvar' object,
%
%       (A*)_{ij} = (A_{ji})*,
%
% following Sec. 8.2 of the sopvar document. The block grid transposes and
% each block is adjoined by @sopvar/ctranspose or @sdopvar/ctranspose.
%
% INPUTS
% - A:      an M x N 'mdopvar' object;
%
% OUTPUTS
% - At:     the N x M adjoint, mapping the output spaces of A to its input
%           spaces;
%
% NOTES
% The container metadata transposes with the grid: the output spaces of At
% are the input spaces of A and vice versa, so space_out and space_in swap,
% as do dim_out and dim_in. The variable registry, the domains and the
% decision variable list are unchanged - an adjoint neither introduces a
% spatial variable nor a decision variable.
%
% Zero blocks stay zero: the adjoint of an absent block is absent, and the
% container metadata continues to determine what it maps between.
%
% The per-block adjoint is where the subtlety lives, and it is already
% handled by the block classes: in an integral direction ZL and ZR swap, but
% in a multiplier direction they must NOT, because delta(s_k - s_k')
% collapses the kernel to the diagonal. See the CANONICAL MULTIPLIER FORM
% note in 'sopvar'. Nothing at the container level may assume the two bases
% exchange.
%
% Cost: M*N block adjoints. No basis is unioned and no decision variable
% list is touched, so this is the cheapest of the three operations. Sec.
% 8.3.1 observes that a row factorization of A becomes a column
% factorization of A*, which is correct, and is one reason those
% factorizations belong inside the operation that wants them rather than in
% storage.
%
% See also PLUS, MTIMES, MDOPVAR, MOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ctranspose(mdopvar)
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
% Initial coding MMP, 09/17/2026. Split out of @mopvar; identical except for
%                the container class it builds and the Zd it carries through.

if ~isa(A,'mdopvar')
    error('ctranspose:badInput','Input must be a mdopvar object.')
end

[M,N] = size(A);
CA = A.C;
Ct = cell(N,M);
for i = 1:M
    for j = 1:N
        if isempty(CA{i,j}),    continue,   end
        Ct{j,i} = CA{i,j}';
    end
end

meta = metadata(A);
meta.space_out = A.space_in;        meta.space_in = A.space_out;
meta.dim_out = A.dim_in;            meta.dim_in = A.dim_out;
At = mdopvar(Ct,meta);

end
