function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = PLUS(A,B) returns the sum A + B of two 'sopvar' objects that map
% between the same spaces, on the unions of their monomial bases.
%
% INPUT
% - A, B:   'sopvar' objects with the same dims, spatial variables and
%           domains. One of them may be a 'dpvar' object or a numeric
%           scalar: a scalar c stands for c*I on a square block, a matrix
%           of the size of the block for its multiplier;
%
% OUTPUT
% - C:      'sopvar' object representing A + B, with ZL and ZR the unions
%           of the monomial bases of A and B; an 'sdopvar' object when a
%           summand is a 'dpvar' with decision variables.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - plus
%
% Copyright (C)2026  M. Peet, S. Shivakumar
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
% Initial coding MMP, SS  - 1_16_2026

% A 'dpvar' or numeric scalar summand, before the checks below: a scalar c
% is c*I, a matrix the multiplier of the block.
if isa(A,'dpvar') || isa(B,'dpvar') || ...
        (isnumeric(A) && isscalar(A)) || (isnumeric(B) && isscalar(B))
    C = dpvar_op_copvar('plus',A,B);
    return
end

% The summands must have the same dimensions, variables, domains and cells
if any(A.dims~=B.dims)
    error('Summands A and B have different dimensions');
end

% The same spatial variables, in the same order, on both sides; lists of
% different length differ
if ~isequal(A.vars.in(:),B.vars.in(:)) || ~isequal(A.vars.out(:),B.vars.out(:))
    error('Summands A and B map between different spaces');
end

if any(any(A.dom.in~=B.dom.in)) || any(any(A.dom.out~=B.dom.out))
    error('input or output variables in summands A and B have different domains');
end

if numel(A.params)~=numel(B.params)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end

% Union monomial bases Z of the two summands, with Zi = Ci*Z the embedding
% of the basis of each summand
ZL1 = A.ZL;
ZL2 = B.ZL;
ZR1 = A.ZR;
ZR2 = B.ZR;

[ZR, C1R, C2R] = UnionBasisMonomials(ZR1, ZR2);
[ZL, C1L, C2L] = UnionBasisMonomials(ZL1, ZL2);

% Lift the embeddings over the matrix dimensions: rows and columns of a
% parameter are component outer, monomial inner
C1L = kron(eye(A.dims(1)), C1L);
C2L = kron(eye(B.dims(1)), C2L);
C1R = kron(eye(A.dims(2)), C1R);
C2R = kron(eye(B.dims(2)), C2R);

% [] and scalar 0 are zero blocks of the size of the block
for i=1:numel(A.params)
    if isempty(A.params{i}) || (isscalar(A.params{i}) && A.params{i}==0)
        A.params{i} = sparse(size(C1L,1),size(C1R,1));
    end
    if isempty(B.params{i}) || (isscalar(B.params{i}) && B.params{i}==0)
        B.params{i} = sparse(size(C2L,1),size(C2R,1));
    end
end

% Add the parameters cell by cell on the union bases (linear indexing of
% the 3 x ... x 3 cell array)
params = A.params;
for i=1:numel(A.params)
    params{i} = C1L'*A.params{i}*C1R +C2L'*B.params{i}*C2R;
end

% Declare the operator (the constructor enforces the canonical form)
C = sopvar(params, A.vars, ZL, ZR, A.dom, A.dims);

end
