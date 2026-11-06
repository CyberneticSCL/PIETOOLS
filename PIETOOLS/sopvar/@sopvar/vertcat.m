function [Pcat] = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCAT = VERTCAT(A_1,...,A_n) returns the vertical concatenation
% [A_1;...;A_n] of 'sopvar' objects on the same spaces: the operator whose
% i-th block of output components is A_i applied to the whole input.
%
% INPUT
% - A_1,...,A_n: two or more 'sopvar' objects with the same number of
%               columns (dims(2)), spatial variables and domains. A 'dpvar'
%               entry is also accepted, in which case the call is routed to
%               'dpvar_op_copvar';
%
% OUTPUT
% - Pcat:       'sopvar' object representing [A_1;...;A_n], with dims(1)
%               the sum of those of the inputs, and with ZL, ZR the unions
%               of the monomial bases of the inputs. A container if any
%               input is a 'dpvar'.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
% or S. Shivakumar at sshivak8@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETools - vertcat
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
% MMP, 09/29/2026: Three or more operands recursed with horzcat, so [A;B;C]
%                  returned [[A;B],C] (a wrong operator, silently, when the
%                  shapes allowed it). Recurse with vertcat.
% MMP, 09/29/2026: A dpvar operand is routed to
%                  'dpvar_op_copvar' (legacy @opvar/vertcat semantics: the
%                  entry maps the blocks' input space into R^{rows}, an
%                  integral from L2), which returns a CONTAINER, since that
%                  entry's output space is not this block's. Before, it
%                  errored "Currently supported only for sopvar".
% MMP, 11/06/2026: Comment pass (maintainer, 11/06/2026): help block in the
%                  house format (the criteria list referred to the fields
%                  of another class), step comments, stale comments removed
%                  (an orphaned 'matrix' remark in the class check, a dead
%                  'Pcat.dim' line); code unchanged.

% A single operand is returned as is
if nargin==1
    Pcat = varargin{1};
    return
end
% dpvar entry: an input-space -> R^{rows} block, so the result              % MMP, 09/29/2026
% is a container; see 'dpvar_op_copvar'.                                    % MMP, 09/29/2026
if any(cellfun(@(a) isa(a,'dpvar'),varargin))                               % MMP, 09/29/2026
    Pcat = dpvar_op_copvar('vertcat',varargin{:});                          % MMP, 09/29/2026
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

% The input dimensions must match
if any(a.dims(2)~=b.dims(2))
    error('Cannot concatenate vertically: Input dimensions of sopvar objects do not match')
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

% Each parameter on the union bases is the vertical concatenation of the
% two re-expressed parameters
params = a.params;
for ii=1:numel(a.params)
    params{ii} = [C1L'*a.params{ii}*C1R; C2L'*b.params{ii}*C2R];
end

dims = a.dims;
dims(1) = dims(1) + b.dims(1);

% Declare the operator (the constructor enforces the canonical form)
Pcat = sopvar(params, a.vars, ZL, ZR, a.dom, dims);                         % MMP, 08/29/2026

% Remaining operands
if nargin>2
%   Pcat = horzcat(Pcat, varargin{3:end});                                  % MMP, 09/29/2026 (was)
    Pcat = vertcat(Pcat, varargin{3:end});                                  % MMP, 09/29/2026
end


end
