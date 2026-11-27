function [Pcat] = horzcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCAT = HORZCAT(A_1,...,A_n) returns the horizontal concatenation
% [A_1,...,A_n] of 'sopvar' objects on the same spaces: the operator that
% applies A_i to the i-th block of input components and sums the results.
%
% INPUT
% - A_1,...,A_n: two or more 'sopvar' objects with the same number of rows
%               (dims(1)), spatial variables and domains. A 'dpvar' entry
%               is also accepted; it maps R^{cols} into the output space;
%
% OUTPUT
% - Pcat:       'sopvar' object representing [A_1,...,A_n], with dims(2)
%               the sum of those of the inputs, and with ZL, ZR the unions
%               of the monomial bases of the inputs. A container if any
%               input is a 'dpvar': 'dopvar' when the 'dpvar' has
%               decision variables, 'opvar' otherwise.
%
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

% A single operand is returned as is
if nargin==1
    Pcat = varargin{1};
    return
end
% dpvar entry: an R^{cols} -> output-space block, so the result
% is a container.
if any(cellfun(@(a) isa(a,'dpvar'),varargin))
    Pcat = dpvar_op_opvar('horzcat',varargin{:});
    return
end

% The first two operands; the rest are handled by recursion below
a = varargin{1};    b = varargin{2};

% Only 'sopvar' operands are supported
if ~isa(b,'sopvar') || ~isa(a, 'sopvar')
    error("Currently supported only for sopvar");
end

% The domains must match
if any(any(a.dom.in~=b.dom.in)) || any(any(a.dom.out~=b.dom.out))
    error('Operators being concatenated have different intervals');
end

% The variables must match. 'isequal', not 'any(~strcmp(...))': for lists
% of different length strcmp gives an empty logical, so {} against {'s1'}
% would pass and an L2 component would be retyped as R.
if ~isequal(a.vars.in(:),b.vars.in(:)) || ~isequal(a.vars.out(:),b.vars.out(:))
    error('Operators being concatenated have input/output variables');
end

% The output dimensions must match
if any(a.dims(1)~=b.dims(1))
    error('Cannot concatenate horizontally: Output dimensions of sopvar objects do not match')
end

% Union monomial bases Z of the two operands, with Zi = Ci*Z the embedding
% of the basis of each operand
ZL1 = a.ZL;
ZL2 = b.ZL;
ZR1 = a.ZR;
ZR2 = b.ZR;

[ZR, C1R, C2R] = UnionBasisMonomials(ZR1, ZR2);
[ZL, C1L, C2L] = UnionBasisMonomials(ZL1, ZL2);

% Lift the embeddings over the matrix dimensions: rows and columns of a
% parameter are component outer, monomial inner
C1L = kron(eye(a.dims(1)), C1L);
C2L = kron(eye(b.dims(1)), C2L);
C1R = kron(eye(a.dims(2)), C1R);
C2R = kron(eye(b.dims(2)), C2R);

% Each parameter on the union bases is the horizontal concatenation of the
% two re-expressed parameters
params = a.params;
for ii=1:numel(a.params)
    params{ii} = [C1L'*a.params{ii}*C1R C2L'*b.params{ii}*C2R];
end

dims = a.dims;
dims(2) = dims(2) + b.dims(2);

% Declare the operator (the constructor enforces the canonical form)
Pcat = sopvar(params, a.vars, ZL, ZR, a.dom, dims);

% Remaining operands
if nargin>2
    Pcat = horzcat(Pcat, varargin{3:end});
end


end
