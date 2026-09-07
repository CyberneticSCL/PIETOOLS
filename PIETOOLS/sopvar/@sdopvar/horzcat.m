function Pcat = horzcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pcat = horzcat(A,B,...) or [A,B,...] concatenates sdopvar operators
% horizontally, so that Pcat maps the stacked input (x_A;x_B;...) to
% A*x_A + B*x_B + ... . All operands must share the output space.
%
% INPUT
% varargin: sdopvar objects with matching vars, dom and output dimension. A
%           fixed 'sopvar' operand is accepted and promoted;
%
% OUTPUT
% Pcat:     sdopvar object of dimension [A.dims(1), sum of input dimensions]
%
% NOTES:
% All N operands are synchronized in ONE pass by 'sync_basis' and then
% concatenated once. Recursing pairwise instead would re-map operand 1
% through N-1 successive transforms, each acting on the growing accumulated
% operator; measured that way the total cost grew as about k^1.56 in the
% number of operands rather than linearly.
%
% Once the operands share bases the concatenation itself is free: columns of
% the coefficient matrix C_gamma are ordered with the matrix column index
% outer and the ZR monomial index inner, and vec is column-major, so
% appending input components appends to vec(C_gamma). Hence the coefficients
% stack along A and sit side by side in B, whose columns are the vec
% positions. The decision variables are the rows of B and are untouched.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
% Initial coding MMP, 09/06/2026
% MMP, 09/07/2026: Accept a fixed 'sopvar' operand, promoting it with
%                  'sopvar2sdopvar'. Building a decision operator out of a
%                  mix of fixed and decision blocks is the normal usage.
% MMP, 09/07/2026: Synchronize all N operands in one pass via 'sync_basis'
%                  instead of recursing pairwise, which re-synchronized the
%                  growing accumulator once per operand.

% Deal with single input case
if nargin==1
    Pcat = varargin{1};
    return
end

% A fixed 'sopvar' operand is promoted to a decision operator with a zero B.
% The decision variable list is taken from whichever operand already has
% one, so 'sync_basis' below sees a single common list.
iss = cellfun(@(x) isa(x,'sopvar'),varargin);
if any(iss)
    isd = find(cellfun(@(x) isa(x,'sdopvar'),varargin),1);
    if isempty(isd),    Zd_p = cell(0,1);
    else,               Zd_p = varargin{isd}.Zd;
    end
    for ii = find(iss)
        varargin{ii} = sopvar2sdopvar(varargin{ii},Zd_p);
    end
end
ops = varargin;
N = numel(ops);

% % % Error handling: every operand must be compatible with the first
if ~all(cellfun(@(x) isa(x,'sdopvar'),ops))
    error('Concatenation is currently supported only between sdopvar and sopvar objects');
end
for k = 2:N
    % Compared with 'isequal' rather than 'strcmp': strcmp of two cellstr of
    % different length returns an empty array, so 'any(~strcmp(...))' is
    % false and a mismatch between {} and {'s1'} would pass unnoticed.
    if ~isequal(ops{k}.vars.in(:),ops{1}.vars.in(:)) || ...
       ~isequal(ops{k}.vars.out(:),ops{1}.vars.out(:))
        error('Operators being concatenated map between different spaces');
    end
    if any(any(ops{k}.dom.in~=ops{1}.dom.in)) || any(any(ops{k}.dom.out~=ops{1}.dom.out))
        error('Operators being concatenated have different intervals');
    end
    if ops{k}.dims(1)~=ops{1}.dims(1)
        error('Cannot concatenate horizontally: output dimensions of sdopvar objects do not match');
    end
    if numel(ops{k}.params.A)~=numel(ops{1}.params.A)
        error('number of terms in the operands is not equal -- one of them is probably malformed');
    end
end

% % % One synchronization for all N operands
[ops,T,Zd,ZL,ZR] = sync_basis(ops);

% % % Append the input components. Each operand contributes a contiguous
% run of vec entries, so A stacks and B sits side by side.
ncell = numel(ops{1}.params.A);
params_new.A = cell(size(ops{1}.params.A));
params_new.B = cell(size(ops{1}.params.B));
Ak = cell(1,N);     Bk = cell(1,N);
for i = 1:ncell
    for k = 1:N
        [Ak{k},Bk{k}] = apply_basis_map(T{k},ops{k}.params.A{i},ops{k}.params.B{i});
    end
    % 'cat' rather than the bracket form, so this cannot re-enter the method
    params_new.A{i} = cat(1,Ak{:});
    params_new.B{i} = cat(2,Bk{:});
end

dims = [ops{1}.dims(1), sum(cellfun(@(P) P.dims(2),ops))];
Pcat = sdopvar(params_new,ops{1}.vars,Zd,ZL,ZR,ops{1}.dom,dims);

end
