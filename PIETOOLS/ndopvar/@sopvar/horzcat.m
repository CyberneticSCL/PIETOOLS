function [Pcat] = horzcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Pcat] = horzcat(varargin) takes n inputs and concatentates them vertically,
% provided they satisfy the following criteria:
% 1) All input must be of type 'sopvar' 
 
% NOTES:
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - horzcat
%
% Copyright (C)2019  M. Peet, S. Shivakumar
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial AT 01/21/2026
% Update to new sopvar AT 05/18/26

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% Extract the operators
a = varargin{1};    b = varargin{2};

% Currently support some matrix-opvar concatenation
if ~isa(b,'sopvar') || ~isa(a, 'sopvar')
    error("Currently supported only for sopvar");
    % Only supported if a in fact corresponds to a matrix as well
end


% Check that domain and variables match
if any(any(a.dom.in~=b.dom.in)) || any(any(a.dom.out~=b.dom.out))
    error('Operators being concatenated have different intervals');
end

% Check that domain and variables match
if any(~strcmp(a.vars.in, b.vars.in))|| any(~strcmp(a.vars.out, b.vars.out))
    error('Operators being concatenated have input/output variables');
end
 

% Check that the output dimensions match
if any(a.dims(1)~=b.dims(1))
    error('Cannot concatenate horizontally: Output dimensions of sopvar objects do not match')
end  

ZL1 = a.ZL;
ZL2 = b.ZL;
ZR1 = a.ZR;
ZR2 = b.ZR;

[ZR, C1R, C2R] = UnionBasisMonomials(ZR1, ZR2);
[ZL, C1L, C2L] = UnionBasisMonomials(ZL1, ZL2);


C1L = kron(eye(a.dims(1)), C1L);
C2L = kron(eye(b.dims(1)), C2L);
C1R = kron(eye(a.dims(2)), C1R);
C2R = kron(eye(b.dims(2)), C2R);

% Finally, let's actually concatenate
params = a.params;
for ii=1:numel(a.params)
    params{ii} = [C1L'*a.params{ii}*C1R C2L'*b.params{ii}*C2R]; % concatenation of quadpoly
end

% Pcat.dim = [a.dim(1), a.dim(2)+b.dim(2)];


dims = a.dims;
dims(2) = dims(2) + b.dims(2);

Pcat = sopvar(params, a.vars, ZL, ZR, a.dom, dims);


if nargin>2 
    Pcat = horzcat(Pcat, varargin{3:end});
end


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