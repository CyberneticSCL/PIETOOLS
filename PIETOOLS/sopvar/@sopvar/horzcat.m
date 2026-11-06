function [Pcat] = horzcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCAT = HORZCAT(A_1,...,A_n) returns the horizontal concatenation
% [A_1,...,A_n] of 'sopvar' objects on the same spaces: the operator that
% applies A_i to the i-th block of input components and sums the results.
%
% INPUT
% - A_1,...,A_n: two or more 'sopvar' objects with the same number of rows
%               (dims(1)), spatial variables and domains. A 'dpvar' entry
%               is also accepted, in which case the call is routed to
%               'dpvar_op_copvar';
%
% OUTPUT
% - Pcat:       'sopvar' object representing [A_1,...,A_n], with dims(2)
%               the sum of those of the inputs, and with ZL, ZR the unions
%               of the monomial bases of the inputs. A container if any
%               input is a 'dpvar'.
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
% Update to new sopvar AT 05/18/26
% Correct order ZL and ZR, split off UnionBasisMonomials, DJ 06/08/2026
% MMP, 09/07/2026: Compare the spatial variable lists with 'isequal'
%                  instead of 'any(~strcmp(...))'. strcmp of two cellstr of
%                  different length scalar-expands to an empty logical, so
%                  the guard was false for a genuine mismatch and blocks on
%                  different spaces were accepted, discarding a spatial
%                  variable and retyping an L2 component as R.
% MMP, 09/29/2026: A dpvar operand is routed to
%                  'dpvar_op_copvar' (legacy @opvar/horzcat semantics: the
%                  entry maps R^{cols} into the blocks' output space), which
%                  returns a CONTAINER, since that entry's input space is not
%                  this block's. Before, it errored "Currently supported only
%                  for sopvar".
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): help block in the
%                  house format, step comments, stale comments removed (a
%                  dangling criteria line in the help text, an orphaned
%                  'matrix' remark in the class check, a dead 'Pcat.dim'
%                  line); code unchanged.

% A single operand is returned as is
if nargin==1
    Pcat = varargin{1};
    return
end
% dpvar entry: an R^{cols} -> output-space block, so the result             % MMP, 09/29/2026
% is a container; see 'dpvar_op_copvar'.                                    % MMP, 09/29/2026
if any(cellfun(@(a) isa(a,'dpvar'),varargin))                               % MMP, 09/29/2026
    Pcat = dpvar_op_copvar('horzcat',varargin{:});                          % MMP, 09/29/2026
    return                                                                  % MMP, 09/29/2026
end                                                                         % MMP, 09/29/2026

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

% The variables must match. 'isequal', not 'any(~strcmp(...))': for lists   % MMP, 09/07/2026
% of different length strcmp gives an empty logical, so {} against {'s1'}   % MMP, 09/07/2026
% passed and an L2 component was silently retyped as R.                     % MMP, 09/07/2026
if ~isequal(a.vars.in(:),b.vars.in(:)) || ~isequal(a.vars.out(:),b.vars.out(:)) % MMP, 09/07/2026
%if any(~strcmp(a.vars.in, b.vars.in))|| any(~strcmp(a.vars.out, b.vars.out)) % MMP, 09/07/2026 (was)
    error('Operators being concatenated have input/output variables');
end

% The output dimensions must match
if any(a.dims(1)~=b.dims(1))
    error('Cannot concatenate horizontally: Output dimensions of sopvar objects do not match')
end

% Union monomial bases Z of the two operands, with Zi = Ci*Z the embedding
% of each operand's basis
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
Pcat = sopvar(params, a.vars, ZL, ZR, a.dom, dims);                         % MMP, 08/29/2026

% Remaining operands
if nargin>2
    Pcat = horzcat(Pcat, varargin{3:end});
end


end
