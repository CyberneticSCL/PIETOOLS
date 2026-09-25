function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = MTIMES(A,B) composes two 'copvar' objects,
%
%       (A*B)_{ij} = sum_{k=1}^{K} A_{ik} * B_{kj},
%
% following Sec. 8.2 of the sopvar document. Also handles a scalar factor.
%
% INPUTS
% - A:      M x K 'copvar', or a scalar double;
% - B:      K x N 'copvar', or a scalar double; the input spaces of A must
%           be the output spaces of B;
% OUTPUTS
% - C:      the M x N composition;
%
% NOTES
% Sec. 8.3.1 writes the inner sum as running to N and "each consists of N
% compositions"; both are K, the shared block dimension, as Sec. 8.2 has it.
% The stated total K*M*N is right.
%
% Absent blocks are skipped rather than composed against, so a structurally
% sparse operator - which a PIE is - costs only its populated products. A
% product that receives no term stays [] and the container metadata says
% what that zero block maps between.
%
% The inner sum folds with 'plus'. There is nothing to batch: 'sopvar'
% carries no decision variables, so the synchronization that makes
% 'plus_batch' worthwhile for a decision container does not arise here. See
% @cdopvar/mtimes, which does batch.
%
% 'copvar' is closed under composition. A decision container is not, which
% is why @cdopvar/mtimes refuses two decision factors.
%
% COST - READ BEFORE SCALING UP. K*M*N block products, and the per-product
% cost is NOT constant: each runs the Sec. 4 composition, whose alpha,beta
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
% block algebra, not in this loop: Sec. 8.3.1's row/column factorization
% removes the repeated basis alignment (K*M*N -> M*K + K*N) but leaves the
% K*M*N kernel-algebra calls, at the inflated union inner dimension. If it
% is added it belongs here as a transient, not in storage.
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
% MMP, 09/25/2026: Renamed the container classes mopvar -> copvar and
%                  mdopvar -> cdopvar, with every file and function named after
%                  them. Mechanical rename, no functional change. Moved from
%                  @mopvar/ with the class.

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

if ~isa(A,'copvar') || ~isa(B,'copvar')
    error('mtimes:badInput','Factors must be copvar objects, or one a scalar.')
end
[M,K] = size(A);        [KB,N] = size(B);
if K~=KB
    error('mtimes:gridMismatch','Inner block dimensions differ: A is %dx%d, B is %dx%d.',...
        M,K,KB,N)
end
if ~isequal(A.vars,B.vars) || ~isequal(A.dom,B.dom)
    error('mtimes:registryMismatch',['Factors are built on different variable '...
        'registries; rebuild them over a common set of variables and domains.'])
end
% The image of B must be the domain of A: input space k of A is output space
% k of B. Compared as masks over the shared registry.
if ~isequal(A.space_in,B.space_out)
    error('mtimes:spaceMismatch','The input spaces of A are not the output spaces of B.')
end
if ~isequal(A.dim_in(:),B.dim_out(:))
    error('mtimes:dimMismatch',...
        'The input dimensions of A do not match the output dimensions of B.')
end
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
