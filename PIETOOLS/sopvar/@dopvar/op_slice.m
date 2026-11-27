function Psop = op_slice(P,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop = op_slice(P,indr,indc) extracts a slice from a decision PI operator
% P: L2^p1[s^1] x ... x L2^pN[s^N] to L2^q1[t^1] x ... x L2^qM[t^M]
%
% INPUT
% P:    dopvar class object to slice
% indr: row indices of the slice: ':', logical, or integers in
%       1:sum(P.dim_out)
% indc: column indices of the slice: ':', logical, or integers in
%       1:sum(P.dim_in)
%
% OUTPUT
% Psop: dopvar object representing P(indr,indc)
%
% NOTES:
% - The indices count rows and columns of the full operator, not blocks:
%   row r of block row i is row sum(P.dim_out(1:i-1))+r, as in size(P);
%   size(P,'block') returns the number of blocks.
% - Each run of consecutive indices within one block row (column) of P
%   becomes one block row (column) of Psop. The order of the indices is
%   kept, and block rows (columns) of P that are not selected are dropped.
% - A block row or column of Psop with only zero blocks gets an explicit
%   zero 'sopvar' block, as a valid container requires one.
% - Psop keeps the decision variables Zd of P, including decision variables
%   that no longer appear in the slice.
% - Psop = op_slice(P,ind) takes linear indices ind, for P with one row or
%   one column. If indr and indc select every row and column in order,
%   Psop = P.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - op_slice (dopvar)
%
% Copyright (C) 2026 PIETOOLS Team
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
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu
%
% DJ, 11/14/2026: Initial coding;

% Numbers of rows and columns of the full operator
nr = sum(P.dim_out);    nc = sum(P.dim_in);
% Linear indices: only a container with one input (output) component in
% all has a sub-operator for them; P(:) of a row would be a reshape.
if isscalar(varargin)
    if nc==1
        varargin = {varargin{1},1};
    elseif nr==1 && ~ischar(varargin{1})
        varargin = {1,varargin{1}};
    else
        error('dopvar:linearIndex',['Linear indexing P(I) needs a dopvar with '...
            'one input component, or one output component and I not '':''; '...
            'use P(I,J).'])
    end
elseif numel(varargin)~=2
    error('dopvar:numIndices','Index a dopvar as P(I,J), or P(I) for one row or column.')
end
% Row and column indices as columns of integers
indr = parse_index(varargin{1},nr,'Row');
indc = parse_index(varargin{2},nc,'Column');

% Every component in order: P itself, no block copied.
if isequal(indr,(1:nr)') && isequal(indc,(1:nc)')
    Psop = P;
    return
end

% Split the row (k=1) and column (k=2) indices into runs of consecutive
% entries in one block of P, each run a block row (column) of the slice:
% blk{k} the block of P of each run, loc{k} its indices inside that block.
ind = {indr,indc};                  d = {P.dim_out,P.dim_in};
blk = {zeros(0,1),zeros(0,1)};      loc = {cell(0,1),cell(0,1)};
for k = 1:2
    if isempty(ind{k})
        continue
    end
    edges = cumsum([0;d{k}(:)]);
    % Block of each index: one more than the number of block ends below it
    % (a block with no components ends where its predecessor does).
    b = sum(ind{k} > edges(2:end-1).',2) + 1;
    first = [true; b(2:end)~=b(1:end-1)];
    blk{k} = b(first);
    loc{k} = mat2cell(ind{k}-edges(b),diff([find(first);numel(b)+1]),1);
end
% Metadata of the slice: the spaces of the selected blocks, one component
% per selected index
[rb,cb] = blk{:};       [rloc,cloc] = loc{:};
meta = metadata(P);
meta.space_out = P.space_out(rb,:);     meta.space_in = P.space_in(cb,:);
meta.dim_out = cellfun(@numel,rloc);    meta.dim_in = cellfun(@numel,cloc);

% Slice each nonzero block (sopvar or sdopvar); a zero block stays [].
C0 = P.C;
C = cell(numel(rb),numel(cb));
for t = 1:numel(rb)
    for u = 1:numel(cb)
        if ~isempty(C0{rb(t),cb(u)})
            C{t,u} = op_slice(C0{rb(t),cb(u)},rloc{t},cloc{u});
        end
    end
end
% A block row or column with no block gets an explicit zero block, at (t,1)
% and then at (1,u).
if ~isempty(C)
    occ = ~cellfun(@isempty,C);
    er = find(~any(occ,2));         occ(er,1) = true;
    ec = find(~any(occ,1))';
    tu = [er,ones(size(er)); ones(size(ec)),ec];
    for k = 1:size(tu,1)
        m1 = meta;
        m1.space_out = meta.space_out(tu(k,1),:);   m1.dim_out = meta.dim_out(tu(k,1));
        m1.space_in = meta.space_in(tu(k,2),:);     m1.dim_in = meta.dim_in(tu(k,2));
        Z = mat2opvar_grid(sparse(m1.dim_out,m1.dim_in),m1,false);
        C{tu(k,1),tu(k,2)} = Z.C{1,1};
    end
end
Psop = dopvar(C,meta);

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function ind = parse_index(ind,n,side)
% ind = parse_index(ind,n,side) returns the subscript ind into 1:n as a
% column of indices.
%
% INPUT
% ind:  ':', logical, or numeric indices
% n:    number of components
% side: 'Row' or 'Column', for the error message
%
% OUTPUT
% ind:  numel x 1 array of integers in 1:n

% ':' selects every component; a logical index becomes integer indices
if ischar(ind) && isequal(ind,':')
    ind = (1:n)';
    return
end
if islogical(ind)
    ind = find(ind);
end
% 'all' of the conditions, so an empty index passes and NaN fails; a
% fractional index would fall between two components.
if ~isnumeric(ind) || ~all(ind(:)>=1 & ind(:)<=n & mod(ind(:),1)==0)
    error('dopvar:badIndex',['%s indices into a dopvar must be '':'', '...
        'logical, or integers in 1..%d (components, not blocks).'],side,n)
end
ind = double(ind(:));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
