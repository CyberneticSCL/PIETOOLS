function At = ctranspose(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% At = CTRANSPOSE(A) returns the adjoint of a 'dopvar' object,
%
%       (A*)_{ij} = (A_{ji})*.
%
% The block grid transposes and each block is replaced by its adjoint.
%
% INPUTS
% - A:      an M x N 'dopvar' object;
%
% OUTPUTS
% - At:     the N x M adjoint, mapping the output spaces of A to its input
%           spaces;
%
% NOTES
% The container metadata transposes with the grid: the output spaces of At
% are the input spaces of A and vice versa, so space_out and space_in swap,
% as do dim_out and dim_in. The variable registry, the domains and the
% decision variable list are unchanged: an adjoint neither introduces a
% spatial variable nor a decision variable.
%
% Zero blocks stay zero: the adjoint of an absent block is absent, and the
% container metadata continues to determine what it maps between.
%
% The per-block adjoint is where the subtlety lives, and it is already
% handled by the block classes: in an integral direction ZL and ZR swap, but
% in a multiplier direction they must NOT, because delta(s_k - s_k')
% collapses the kernel to the diagonal. Nothing at the container level may
% assume the two bases exchange.
%
% Cost: M*N block adjoints. No basis is unioned and no decision variable
% list is touched.
%
% See also PLUS, MTIMES, DOPVAR, OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ctranspose(dopvar)
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
% Initial coding MMP, 09/17/2026

if ~isa(A,'dopvar')
    error('ctranspose:badInput','Input must be a dopvar object.')
end

% Transpose the grid, replacing each nonzero block by its adjoint
[M,N] = size(A,'block');
CA = A.C;
Ct = cell(N,M);
for i = 1:M
    for j = 1:N
        if isempty(CA{i,j}),    continue,   end
        Ct{j,i} = CA{i,j}';
    end
end

% Swap the output and input metadata
meta = metadata(A);
meta.space_out = A.space_in;        meta.space_in = A.space_out;
meta.dim_out = A.dim_in;            meta.dim_in = A.dim_out;
At = dopvar(Ct,meta);

end
