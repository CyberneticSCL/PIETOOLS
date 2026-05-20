function C = plus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = plus(A,B) adds two sopvar operators A,B: L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Date: 1/16/26
% Version: 1.0
% 
% INPUT
% A: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% B: sopvar class object L_2^p[Si,...Sn] to L_2^q[Sj,...,Sm]
% Note: A,B must have identical dimensions and domains
%
% OUTPUT
% C = A+B:  L_2^q[Sj,...,Sm] to L_2^p[Si,...Sn]
% 
% NOTES:
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
% Allowing different monomials AT - 05/26/26

% initialize added output class

% Error handling: Checks to ensure A and B are compatible
if any(A.dims~=B.dims)
    error('Summands A and B have different dimensions');
end

if any(~strcmp(A.vars.in,B.vars.in)) ||any(~strcmp(A.vars.out,B.vars.out))
    error('Summands A and B map between different spaces');
end

if any(any(A.dom.in~=B.dom.in)) || any(any(A.dom.out~=B.dom.out))
    error('input or output variables in summands A and B have different domains');
end

if numel(A.params)~=numel(B.params)
    error('number of terms in summands is not equal -- one of them is probably malformed');
end
% dimensions of params in A and B, vars_in and vars_out


ZL1 = A.ZL;
ZL2 = B.ZL;
ZR1 = A.ZR;
ZR2 = B.ZR;

%find common basis
[ZR, C1R, C2R] = UnionBasisMonomials(ZR1, ZR2);
[ZL, C1L, C2L] = UnionBasisMonomials(ZL1, ZL2);


C1L = kron(eye(A.dims(1)), C1L);
C2L = kron(eye(B.dims(1)), C2L);
C1R = kron(eye(A.dims(2)), C1R);
C2R = kron(eye(B.dims(2)), C2R);



params = A.params;
for i=1:numel(A.params)  % linear indexing of multi-dimensional cell array
    params{i} = C1L'*A.params{i}*C1R +C2L'*B.params{i}*C2R;  % adding quadpoly objects. % Note: chatgpt says this approach is faster than cellfun
end


C = sopvar(params, A.vars, ZL, ZR, A.dom, A.dims);

end




% ---------------- local helpers ----------------
function [Z, C1, C2] = UnionBasisMonomials(Z1, Z2)
% computes monomial basis for to monomial vectors
% and transformation matrix
% such that for any matrix A we have 
% A*Z1 = A*C1*Z and A*Z2 = A*C2*Z
%
% input Z1, Z2 -- cell arrays
% output C1 -- (n1, n_new) sparse matrix
% output C2 -- (n2, n_new) sparse matrix
% output Z -- new basis (cell array)
N = length(Z1);
Z = cell(1:N);

% construct union basis
for i = 1:N
    Z{i} = union(Z1{i}, Z2{i});
end

degmat_1   = Z1{1};
degmat_2   = Z2{1};
degmat_new = Z{1};

% construct full monomial basis
for dim = 2:N
    degmat_new = [kron(degmat_new, ones(length(Z{dim}), 1)),  kron(ones(length(degmat_new), 1), Z{dim})];
    degmat_1 = [kron(degmat_1, ones(length(Z1{dim}), 1)),  kron(ones(length(degmat_1), 1), Z1{dim})];
    degmat_2 = [kron(degmat_2, ones(length(Z2{dim}), 1)),  kron(ones(length(degmat_2), 1), Z2{dim})];
end

% find identical rows 
[~, C1_index] =ismember(degmat_1,  degmat_new, 'rows');
[~, C2_index] =ismember(degmat_2,  degmat_new, 'rows');

n1 = length(degmat_1);
n2 = length(degmat_2);
n_new = length(degmat_new);
C1 = sparse(1:n1, C1_index, ones(n1, 1), n1, n_new);
C2 = sparse(1:n2, C2_index, ones(n2, 1), n2, n_new);
end
