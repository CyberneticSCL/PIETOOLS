function C = mtimes(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = MTIMES(A,B) composes a decision container with a FIXED one,
%
%       (A*B)_{ij} = sum_{k=1}^{K} A_{ik} * B_{kj},
%
% following Sec. 8.2 of the sopvar document. Also handles a scalar factor.
%
% INPUTS
% - A:      M x K 'mdopvar' or 'mopvar', or a scalar double;
% - B:      K x N 'mdopvar' or 'mopvar', or a scalar double; the input
%           spaces of A must be the output spaces of B. AT MOST ONE of A
%           and B may be an 'mdopvar';
% OUTPUTS
% - C:      the M x N composition, an 'mdopvar';
%
% NOTES
% Two 'mdopvar' factors are refused outright, before any work: the product
% of two operators affine in the decision variables is quadratic in them and
% no block class can represent it. This is the check the class split buys -
% it is decided by the argument types, with no scan of the blocks. Compose a
% decision operator only against a fixed one, as in T'*P*T. A container that
% happens to hold no 'sdopvar' block is still refused; use 'mopvar' for the
% fixed factor, which is what it is for.
%
% Sec. 8.3.1 writes the inner sum as running to N and "each consists of N
% compositions"; both are K, the shared block dimension, as Sec. 8.2 has it.
% The stated total K*M*N is right.
%
% Absent blocks are skipped rather than composed against, so a structurally
% sparse operator - which a PIE is - costs only its populated products. A
% product that receives no term stays [] and the container metadata says
% what that zero block maps between.
%
% The inner sum goes to 'plus_batch' when its terms include an 'sdopvar', so
% the K terms are synchronized once rather than K-1 times against a growing
% accumulator; pairwise accumulation profiled at 32% of 'possopvar' at three
% spatial variables. All-'sopvar' terms fold with 'plus' instead, which costs
% nothing extra since they carry no decision variables.
%
% COST - READ BEFORE SCALING UP. K*M*N block products, and the per-product
% cost is NOT constant: each runs the Sec. 4 composition, whose alpha,beta
% double sum is 9^n3 in the number of shared variables. Measured warm on
% 'sopvar' blocks, 2 components each, all variables shared:
%
%   n3        0       1       2       3          4
%   deg 1  0.025 s 0.053 s 0.062 s 0.102 s    1.373 s
%   deg 2  0.025 s 0.025 s 0.023 s 0.139 s   12.537 s
%
% The 9^n3 factor is invisible through n3 = 3, where fixed overhead
% dominates, then takes over. End to end over the subset lattice at degree
% 1, 2 components per space, no decision variables, one composition costs
% 0.86 s at 4 spaces (64 products), 4.97 s at 8 spaces (512) and 56.7 s at
% 16 spaces (4096). With decision variables, a 4x4 grid at q = 2e5 measured
% 5.65 s for its 64 products. This is the scaling wall of the class, and it
% is in the block algebra, not in this loop: Sec. 8.3.1's row/column
% factorization removes the repeated basis alignment (K*M*N -> M*K + K*N)
% but leaves the K*M*N kernel-algebra calls, at the inflated union inner
% dimension. If it is added it belongs here as a transient, not in storage.
%
% See also PLUS, CTRANSPOSE, MDOPVAR, MOPVAR, PLUS_BATCH.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes(mdopvar)
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

% % % Scalar factor: scales every populated block, changes no metadata.
if isnumeric(A) || isnumeric(B)
    if isnumeric(A),    a = A;  P = B;      else,   a = B;  P = A;      end
    if ~isscalar(a)
        error('mtimes:nonScalarNumeric',['A numeric factor must be scalar; '...
            'multiplication by a matrix changes the spaces and must be a container.'])
    end
    Cc = P.C;
    for ii = 1:numel(Cc)
        if ~isempty(Cc{ii}),    Cc{ii} = a*Cc{ii};      end
    end
    C = mdopvar(Cc,metadata(P));
    return
end

% Decided by the argument types, before promotion and before any work: once
% a 'mopvar' operand has been promoted both factors are 'mdopvar', so the
% test has to be made on what the caller actually passed.
if isa(A,'mdopvar') && isa(B,'mdopvar')
    error('mdopvar:decisionTimesDecision',...
        ['Both factors are mdopvar, so the product would be quadratic in the '...
         'decision variables and is not representable. Compose a decision '...
         'operator only against a mopvar, as in T''*P*T.'])
end
if isa(A,'mopvar'),     A = mdopvar(A);     end
if isa(B,'mopvar'),     B = mdopvar(B);     end
if ~isa(A,'mdopvar') || ~isa(B,'mdopvar')
    error('mtimes:badInput',...
        'Factors must be mdopvar or mopvar objects, or one a scalar.')
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

% Reconcile the decision variables once for the whole composition, so that
% none of the K*M*N products merges a name list of its own.
CA = A.C;   CB = B.C;   Zd = A.Zd(:);
if ~isequal(A.Zd(:),B.Zd(:))
    Zd = unique([A.Zd(:);B.Zd(:)],'stable');
    CA = put_on_list(CA,Zd);
    CB = put_on_list(CB,Zd);
end

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
        if nt==1
            Cc{i,j} = terms{1};
        elseif nt>1 && any(cellfun(@(t) isa(t,'sdopvar'),terms(1:nt)))
            Cc{i,j} = plus_batch(terms{1:nt});      % one synchronization
        elseif nt>1
            Cc{i,j} = terms{1};
            for t = 2:nt,   Cc{i,j} = Cc{i,j} + terms{t};    end
        end
    end
end

% The composition maps the input spaces of B to the output spaces of A, both
% already validated, so the metadata is assembled rather than re-derived.
meta = metadata(A);
meta.space_in = B.space_in;     meta.dim_in = B.dim_in;     meta.Zd = Zd;
C = mdopvar(Cc,meta);

end
