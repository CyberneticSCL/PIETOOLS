function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = MTIMES(A,B) composes two 'copvar' objects,
%
%       (A*B)_{ij} = sum_{k=1}^{K} A_{ik} * B_{kj}.
%
% Also handles a scalar factor and a 'dpvar' factor.
%
% INPUTS
% - A:      M x K 'copvar', or a scalar double, or a 'dpvar';
% - B:      K x N 'copvar', or a scalar double, or a 'dpvar'; the input
%           spaces of A must be the output spaces of B. A 'dpvar' scalar
%           scales (gam*P is a 'cdopvar'); a 'dpvar' matrix multiplies the
%           output (D*P) or input (P*D) components. A numeric matrix M is
%           an error; dpvar(M) multiplies by it;
% OUTPUTS
% - C:      the M x N composition;
%
% NOTES
% Absent blocks are skipped rather than composed against, so a structurally
% sparse operator, such as a PIE, costs only its populated products. A
% product that receives no term stays [] and the container metadata says
% what that zero block maps between.
%
% Factors over different registries are first restated over their sorted
% union, so an R^n -> R^n factor such as Dzu, whose registry is empty,
% composes with an L_2 one. Metadata only, O(nv); the result keeps the
% union registry.
%
% The inner sum is added pairwise. There is nothing to batch: 'sopvar'
% carries no decision variables, so no decision variable lists need to be
% synchronized.
%
% 'copvar' is closed under composition. A decision container is not: the
% product of two decision containers is an error.
%
% COST: READ BEFORE SCALING UP. K*M*N block products, and the per-product
% cost is NOT constant: each runs the block composition, whose alpha,beta
% double sum is 9^n3 in the number of shared variables. Measured warm, 2
% components per block, all variables shared:
%
%   n3        0       1       2       3          4
%   deg 1  0.025 s 0.053 s 0.062 s 0.102 s    1.373 s
%   deg 2  0.025 s 0.025 s 0.023 s 0.139 s   12.537 s
%
% The 9^n3 factor is invisible through n3 = 3, where fixed overhead
% dominates, then takes over. End to end over the subset lattice at degree
% 1, 2 components per space, no decision variables, one composition costs
% 0.86 s at 4 spaces (64 products), 4.97 s at 8 spaces (512) and 56.7 s at
% 16 spaces (4096). This is the scaling wall of the class, and it is in the
% block algebra, not in this loop: a row/column factorization would remove
% the repeated basis alignment (K*M*N -> M*K + K*N) but leave the K*M*N
% kernel-algebra calls, at the inflated union inner dimension.
%
% See also PLUS, CTRANSPOSE, COPVAR, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes(copvar)
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

% dpvar factor: a scalar or a matrix multiplier, handled before any
% work, so the composition below is untouched.
if isa(A,'dpvar') || isa(B,'dpvar')
    C = dpvar_op_copvar('mtimes',A,B);
    return
end

% % % Scalar factor: scales every populated block, changes no metadata.
if isnumeric(A) || isnumeric(B)
    if isnumeric(A),    a = A;  P = B;      else,   a = B;  P = A;      end
    if ~isscalar(a)
        error('mtimes:nonScalarNumeric',['A numeric factor must be scalar; '...
            'multiplication by a matrix changes the spaces and must be a copvar.'])
    end
    Cc = P.C;
    for ii = 1:numel(Cc)
        if ~isempty(Cc{ii}),    Cc{ii} = a*Cc{ii};      end
    end
    C = copvar(Cc,metadata(P));
    return
end

% Check the factors and their block grids
if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('mtimes:badInput','Factors must be copvar objects, or one a scalar.')
end
[M,K] = size(A,'block');        [KB,N] = size(B,'block');
if K~=KB
    error('mtimes:gridMismatch','Inner block dimensions differ: A is %dx%d, B is %dx%d.',...
        M,K,KB,N)
end
% One registry for both, so the masks below compare by name; a variable
% on two domains is an error.
[A,B] = merge_copvar_registry('copvar','mtimes',A,B);
% The image of B must be the domain of A: input space k of A is output space
% k of B. Compared as masks over the shared registry.
if ~isequal(A.space_in,B.space_out)
    error('mtimes:spaceMismatch','The input spaces of A are not the output spaces of B.')
end
if ~isequal(A.dim_in(:),B.dim_out(:))
    error('mtimes:dimMismatch',...
        'The input dimensions of A do not match the output dimensions of B.')
end
% Block (i,j) of the product: the sum over k of the nonzero A_{ik}*B_{kj}
CA = A.C;   CB = B.C;
Cc = cell(M,N);
terms = cell(1,K);
for i = 1:M
    for j = 1:N
        nt = 0;
        for k = 1:K
            % A zero block annihilates the product, so skip it. This is the
            % benefit of the zero-block convention.
            if isempty(CA{i,k}) || isempty(CB{k,j}),    continue,   end
            nt = nt+1;      terms{nt} = CA{i,k}*CB{k,j};
        end
        if nt>=1
            Cc{i,j} = terms{1};
            for t = 2:nt,   Cc{i,j} = Cc{i,j} + terms{t};    end
        end
    end
end

% The composition maps the input spaces of B to the output spaces of A, both
% already validated, so the metadata is assembled rather than re-derived.
meta = metadata(A);
meta.space_in = B.space_in;     meta.dim_in = B.dim_in;
C = copvar(Cc,meta);

end
