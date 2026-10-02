function C = plus_batch(varargin)
% C = PLUS_BATCH(P1,P2,...,PN) sums N containers on one block grid,
% reconciling their decision variable lists ONCE, where a chain
% P1 + P2 + ... + PN reconciles once per '+'.
%
% INPUT
% varargin: N >= 1 'cdopvar' / 'copvar' containers, at least one a
%           'cdopvar', on the same block grid, mapping between the same
%           spaces with the same component dimensions (as for '+'). A
%           'copvar' operand is promoted;
%
% OUTPUT
% C:        'cdopvar' equal to P1 + P2 + ... + PN
%
% NOTES
% In a chain, each '+' puts the accumulated sum on the union of its list and
% the next operand's, remapping the rows of every decision block of the sum
% ('setdvars'): the sum's list is a prefix of the union, yet each block is
% rebuilt, so N operands with their own lists cost O(N^2) block rows. Here
% the lists are merged once ('merge_dvar_lists', stable union in operand
% order, the order the chain produces) and each grid cell is summed by one
% sdopvar 'plus_batch' on the shared list, adding the blocks left to right
% as the chain does. A call with 'copvar' operands only dispatches to
% '@copvar/plus_batch'.
% Measured on the heatNd bench face sums, N-D (2N+1 terms), bit-identical to
% the chain: N = 3, q = 2.9e5 and 1.1e6: 4.1 -> 2.0 s, peak commit 1.29 ->
% 1.05 GB; N = 2, q = 1.5e4: 0.25 -> 0.19 s; N = 1: 0.16 -> 0.15 s.
%
% See also PLUS, PLUS_BATCH, MERGE_DVAR_LISTS, MERGE_COPVAR_REGISTRY.
%
% Initial coding MMP, 10/01/2026

ops = varargin;     N = numel(ops);
isc = cellfun(@(x) isa(x,'copvar'),ops);    isd = cellfun(@(x) isa(x,'cdopvar'),ops);
if ~all(isc | isd)
    error('plus_batch:badInput','Summands must be copvar or cdopvar objects.')
end
if N==1,    C = ops{1};     return,     end
for k = find(isc),  ops{k} = cdopvar(ops{k});   end
sz = size(ops{1});
for k = 2:N
    if ~isequal(size(ops{k}),sz)
        error('plus:gridMismatch','Summands have block grids %s and %s.',...
              mat2str(sz),mat2str(size(ops{k})))
    end
end
% One registry for all operands (as '+' does pairwise), then the checks '+'
% makes, against the first operand.
[ops{:}] = merge_copvar_registry('cdopvar','plus',ops{:});
A = ops{1};
for k = 2:N
    if ~isequal(A.space_out,ops{k}.space_out) || ~isequal(A.space_in,ops{k}.space_in)
        error('plus:spaceMismatch','Summands map between different spaces.')
    end
    if ~isequal(A.dim_out(:),ops{k}.dim_out(:)) || ~isequal(A.dim_in(:),ops{k}.dim_in(:))
        error('plus:dimMismatch','Summands have different component dimensions.')
    end
end
% One reconciliation of the decision lists: block cell (ii,k) of operand k
% has source k.
nc = prod(sz);
blocks = cell(1,N*nc);      src = zeros(1,N*nc);    Zds = cell(1,N);
for k = 1:N
    blocks((k-1)*nc+(1:nc)) = reshape(ops{k}.C,1,[]);
    src((k-1)*nc+(1:nc)) = k;
    Zd = ops{k}.Zd;     if ~iscolumn(Zd),   Zd = Zd(:);     end
    Zds{k} = Zd;
end
[blocks,Zd] = merge_dvar_lists(blocks,Zds,src);
% Each grid cell summed once, operands left to right.
Cc = cell(sz);
for ii = 1:nc
    b = blocks(ii:nc:end);              % cell ii of operands 1..N
    b = b(~cellfun(@isempty,b));
    if isempty(b),  continue,   end
    if isscalar(b)
        Cc{ii} = b{1};
    elseif any(cellfun(@(x) isa(x,'sdopvar'),b))
        Cc{ii} = plus_batch(b{:},'shared_Zd');      % every sdopvar is on Zd
    else
        S = b{1};                                   % fixed blocks only
        for t = 2:numel(b),     S = S + b{t};   end
        Cc{ii} = S;
    end
    if isa(Cc{ii},'sdopvar'),   Cc{ii}.Zd = Zd;     end
end
meta = metadata(A);     meta.Zd = Zd;
C = cdopvar(Cc,meta);
end
