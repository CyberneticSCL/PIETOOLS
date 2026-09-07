function Pcat = vertcat(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pcat = vertcat(A,B,...) or [A;B;...] concatenates sdopvar operators
% vertically, so that Pcat maps x to the stacked output (A*x;B*x;...). All
% operands must share the input space.
%
% INPUT
% varargin: sdopvar objects with matching vars, dom and input dimension. A
%           fixed 'sopvar' operand is accepted and promoted;
%
% OUTPUT
% Pcat:     sdopvar object of dimension [sum of output dimensions, A.dims(2)]
%
% NOTES:
% All N operands are synchronized in ONE pass by 'sync_basis' and then
% combined once; see the note in horzcat on why pairwise recursion is worse.
%
% Unlike horzcat this is not a plain append: rows of the coefficient matrix
% C_gamma are ordered with the matrix row index outer and the ZL monomial
% index inner, and vec is column-major, so stacking output components
% interleaves the operands column block by column block. J below holds the
% vec positions operand k occupies, so each contribution is a scatter of
% rows of A and of COLUMNS of B; the decision variables are the rows of B
% and are untouched.
%
% Cross-check: the adjoint of a stacked operator is the concatenation of the
% adjoints, so [A;B] must equal ([A',B'])'.
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
%                  instead of recursing pairwise.
% MMP, 09/07/2026: Scatter each operand with an indexed assignment rather
%                  than collecting all N operands' triplets into one
%                  sparse() call, which was 2x to 4.5x slower because that
%                  constructor sorts every entry.

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
    if ops{k}.dims(2)~=ops{1}.dims(2)
        error('Cannot concatenate vertically: input dimensions of sdopvar objects do not match');
    end
    if numel(ops{k}.params.A)~=numel(ops{1}.params.A)
        error('number of terms in the operands is not equal -- one of them is probably malformed');
    end
end

% % % One synchronization for all N operands
[ops,T,Zd,ZL,ZR] = sync_basis(ops);

% % % Positions each operand occupies in the new vec: within each of the
% ncol columns, operand k takes the rows offset by the operands before it.
NL = prod([cellfun(@numel,ZL),1]);
NR = prod([cellfun(@numel,ZR),1]);
mk   = cellfun(@(P) P.dims(1),ops);
nrow = sum(mk)*NL;      ncol = ops{1}.dims(2)*NR;
off  = [0,cumsum(mk(1:end-1))]*NL;
J = cell(1,N);
for k = 1:N
    J{k} = reshape(off(k)+(1:mk(k)*NL).' + (0:ncol-1)*nrow,[],1);
end

% % % Scatter each operand into its rows. One indexed assignment per        % MMP, 09/07/2026
% operand, not a single triplet sparse() over all N: measured 2x to 4.5x    % MMP, 09/07/2026
% faster across N=2..8 and densities 1e-2..1 (40.9 s against 14.4 s at      % MMP, 09/07/2026
% N=8 with 6e8 nonzeros), because the triplet form must sort every entry.   % MMP, 09/07/2026
q = numel(Zd);
ncell = numel(ops{1}.params.A);
params_new.A = cell(size(ops{1}.params.A));
params_new.B = cell(size(ops{1}.params.B));
for i = 1:ncell
    params_new.A{i} = sparse(nrow*ncol,1);                                  % MMP, 09/07/2026
    params_new.B{i} = sparse(q,nrow*ncol);                                  % MMP, 09/07/2026
    for k = 1:N
        [Ak,Bk] = apply_basis_map(T{k},ops{k}.params.A{i},ops{k}.params.B{i});
        params_new.A{i}(J{k})   = Ak;                                       % MMP, 09/07/2026
        params_new.B{i}(:,J{k}) = Bk;                                       % MMP, 09/07/2026
    end
end
%   ia = cell(1,N);  va = cell(1,N);  ib = cell(1,N);  jb = cell(1,N);      % MMP, 09/07/2026 (was)
%   vb = cell(1,N);                                                         % MMP, 09/07/2026 (was)
%   [r,~,v] = find(Ak);      ia{k} = J{k}(r);   va{k} = v;                  % MMP, 09/07/2026 (was)
%   [r,c,v] = find(Bk);      ib{k} = r; jb{k} = J{k}(c); vb{k} = v;         % MMP, 09/07/2026 (was)
%   params_new.A{i} = sparse(cat(1,ia{:}),1,cat(1,va{:}),nrow*ncol,1);      % MMP, 09/07/2026 (was)
%   params_new.B{i} = sparse(cat(1,ib{:}),cat(1,jb{:}),cat(1,vb{:}),q,nrow*ncol); % MMP, 09/07/2026 (was)

dims = [sum(mk), ops{1}.dims(2)];
Pcat = sdopvar(params_new,ops{1}.vars,Zd,ZL,ZR,ops{1}.dom,dims);

end
