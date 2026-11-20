function C = plus_batch(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% C = plus_batch(P1,P2,...,PN) sums N sdopvar operators, synchronizing all N
% once instead of once per addition.
%
% INPUTS
% varargin: N 'sdopvar' objects of equal dimension, sharing vars and dom. A
%           fixed 'sopvar' operand is accepted and promoted;
%           a trailing 'shared_Zd' is the caller's guarantee that every
%           'sdopvar' operand carries the same decision variable list (e.g.
%           all were built from one Zd); it skips an O(q) comparison per
%           operand. Pass it only when that holds by construction;
%
% OUTPUTS
% C:        'sdopvar' object equal to P1 + P2 + ... + PN
%
% NOTES:
% Equivalent to adding the operands pairwise, but the reconciliation of
% decision variables and monomial bases happens once for all N rather than
% N-1 times against a growing accumulator. Pairwise accumulation was
% measured at 32% of the runtime of a positive-operator declaration in
% three spatial variables, where the 27 basis operators give 729 additions.
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

% A trailing 'shared_Zd' is removed before the operands are processed; see
% INPUT. Any other text input is an error rather than an operand.
shared_Zd = false;
if ~isempty(varargin) && (ischar(varargin{end}) || isstring(varargin{end}))
    if ~strcmp(varargin{end},'shared_Zd')
        error("plus_batch: unrecognized option '%s'.",varargin{end})
    end
    shared_Zd = true;
    varargin(end) = [];
end
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
% A promoted 'sopvar' took Zd_p from an 'sdopvar' operand above, so the
% caller's guarantee covers it too.
[ops,T,Zd,ZL,ZR] = sync_basis(ops,shared_Zd);

% % % Add the summands in sequence. Only the SYNCHRONIZATION is batched;
% the summation stays as sparse additions, which merge two sorted column
% structures (one sparse() over all N triplets must sort N*nnz entries and
% was measured 1.6x slower at N=729, q=2000).
ncell = numel(ops{1}.params.A);
params_new.A = cell(size(ops{1}.params.A));
params_new.B = cell(size(ops{1}.params.B));
for i = 1:ncell
    [Asum,Bsum] = apply_basis_map(T{1},ops{1}.params.A{i},ops{1}.params.B{i});
    for k = 2:N
        % A summand without entries maps to sparse zeros of the merged
        % size, built directly instead of by apply_basis_map's scatter; the
        % guard admits exactly the inputs on which that scatter returns
        % these zeros (real double, numel/columns matching T.idx).
        Ar = ops{k}.params.A{i};    Br = ops{k}.params.B{i};
        if ~isempty(T{k}) && nnz(Ar)==0 && nnz(Br)==0 && isa(Ar,'double') ...
            && isa(Br,'double') && isreal(Ar) && isreal(Br) && ismatrix(Br) ...
            && numel(Ar)==numel(T{k}.idx) && size(Br,2)==numel(T{k}.idx)
            % Skip before building the zeros when 'noop_sum' would skip
            % both adds anyway: running sums real sparse double of the
            % merged sizes (a full sum could hold -0, which + 0 flips).
            if issparse(Asum) && issparse(Bsum) && isa(Asum,'double') ...
                && isa(Bsum,'double') && isreal(Asum) && isreal(Bsum) ...
                && isequal(size(Asum),[T{k}.nC,1]) ...
                && isequal(size(Bsum),[size(Br,1),T{k}.nC])
                continue
            end
            Ak = sparse(T{k}.nC,1);     Bk = sparse(size(Br,1),T{k}.nC);
        else
            [Ak,Bk] = apply_basis_map(T{k},Ar,Br);
        end
        % Adding an empty summand would copy the whole running sum; skip
        % it where the sum keeps the same bits ('noop_sum').
        if ~noop_sum(Asum,Ak),  Asum = Asum + Ak;   end
        if ~noop_sum(Bsum,Bk),  Bsum = Bsum + Bk;   end
    end
    params_new.A{i} = Asum;
    params_new.B{i} = Bsum;
end

C = sdopvar(params_new,ops{1}.vars,Zd,ZL,ZR,ops{1}.dom,ops{1}.dims);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf = noop_sum(S,Z)
% TF = NOOP_SUM(S,Z) is true when S+Z equals S bit for bit, so the add can
% be skipped: both real sparse double of one size, and Z stores nothing.
% Sparse storage holds no zeros, so every value of S passes through as is
% (x+0 == x for x ~= 0); a full S could hold -0, which +0 would flip; a
% class, complexity or size difference would change the result or error.
tf = issparse(S) && issparse(Z) && isa(S,'double') && isa(Z,'double') ...
     && isreal(S) && isreal(Z) && nnz(Z)==0 && isequal(size(S),size(Z));
end
