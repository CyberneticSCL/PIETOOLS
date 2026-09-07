function C = plus_batch(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = plus_batch(P1,P2,...,PN) sums N sdopvar operators, synchronizing all N
% once instead of once per addition.
%
% INPUT
% varargin: N 'sdopvar' objects of equal dimension, sharing vars and dom. A
%           fixed 'sopvar' operand is accepted and promoted;
%
% OUTPUT
% C:        'sdopvar' object equal to P1 + P2 + ... + PN
%
% NOTES:
% Equivalent to folding @sdopvar/plus over the list, but the reconciliation
% of decision variables and monomial bases happens once for all N rather
% than N-1 times against a growing accumulator. Accumulating pairwise was
% measured at 32% of 'possopvar' runtime for three spatial variables, where
% the 27 basis operators give 729 additions.
%
% Only the synchronization is batched. The summation itself is left as
% ordinary sparse additions: batching that too, by concatenating all N
% operands' triplets into a single sparse() call, was measured 1.6x SLOWER
% at N=729, q=2000, because that constructor must sort N*nnz entries while
% sparse addition merges two already-sorted column structures. The decision
% variables are the rows of B and are only ever indexed, never densified.
%
% See also PLUS, SYNC_BASIS.
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
% Initial coding MMP, 09/07/2026

ops = varargin;
if isempty(ops)
    error('plus_batch requires at least one operand.')
elseif isscalar(ops)
    C = ops{1};
    return
end

% A fixed 'sopvar' summand is promoted to a decision operator with a zero B,
% taking its decision variable list from whichever operand already has one.
iss = cellfun(@(x) isa(x,'sopvar'),ops);
if any(iss)
    isd = find(cellfun(@(x) isa(x,'sdopvar'),ops),1);
    if isempty(isd),    Zd_p = cell(0,1);
    else,               Zd_p = ops{isd}.Zd;
    end
    for ii = find(iss)
        ops{ii} = sopvar2sdopvar(ops{ii},Zd_p);
    end
end
N = numel(ops);

% % % Error handling: every summand must be compatible with the first
if ~all(cellfun(@(x) isa(x,'sdopvar'),ops))
    error('plus_batch is supported only for sdopvar and sopvar objects.')
end
for k = 2:N
    if any(ops{k}.dims~=ops{1}.dims)
        error('Dimensions of summands do not match');
    end
    % 'isequal' rather than 'strcmp': strcmp of two cellstr of different
    % length returns an empty array, so any(~strcmp(...)) is false and a
    % mismatch between {} and {'s1'} would pass unnoticed.
    if ~isequal(ops{k}.vars.in(:),ops{1}.vars.in(:)) || ...
       ~isequal(ops{k}.vars.out(:),ops{1}.vars.out(:))
        error('Summands map between different spaces');
    end
    if any(any(ops{k}.dom.in~=ops{1}.dom.in)) || any(any(ops{k}.dom.out~=ops{1}.dom.out))
        error('Input or output variables in the summands have different domains');
    end
    if numel(ops{k}.params.A)~=numel(ops{1}.params.A)
        error('number of terms in summands is not equal -- one of them is probably malformed');
    end
end

% % % One synchronization for all N summands
[ops,T,Zd,ZL,ZR] = sync_basis(ops);

% % % Add the summands in sequence. Only the SYNCHRONIZATION is batched;
% the summation is left as sparse additions, which merge two already-sorted
% column structures. Collecting all N operands' triplets into one sparse()
% call instead was measured 1.6x SLOWER at N=729, q=2000 (16.6 s against
% 10.1 s), because that constructor must sort N*nnz entries; it was faster
% only when q was small enough for the sort to be free.
ncell = numel(ops{1}.params.A);
params_new.A = cell(size(ops{1}.params.A));
params_new.B = cell(size(ops{1}.params.B));
for i = 1:ncell
    [Asum,Bsum] = apply_basis_map(T{1},ops{1}.params.A{i},ops{1}.params.B{i});
    for k = 2:N
        [Ak,Bk] = apply_basis_map(T{k},ops{k}.params.A{i},ops{k}.params.B{i});
        Asum = Asum + Ak;
        Bsum = Bsum + Bk;
    end
    params_new.A{i} = Asum;
    params_new.B{i} = Bsum;
end

C = sdopvar(params_new,ops{1}.vars,Zd,ZL,ZR,ops{1}.dom,ops{1}.dims);

end
