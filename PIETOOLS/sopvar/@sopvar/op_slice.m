function Psop = op_slice(P,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Psop = op_slice(P,indr,indc) extracts a slice from a PI operator
% P: L2^p[S1,S3] to L2^q[S2,S3]
%
% INPUT
% P:    sopvar class object to slice
% indr: row indices of the slice: ':', logical, or integers in 1:q
% indc: column indices of the slice: ':', logical, or integers in 1:p
%
% OUTPUT
% Psop: sopvar object representing P(indr,indc)
%
% NOTES:
% - Psop = op_slice(P,ind) takes linear indices ind, for P with one row or
%   one column. P(:) of a row operator is refused, as it is not a slice.
% - Psop keeps the monomial bases ZL and ZR of P.
% - If indr and indc select every row and column in order, Psop = P.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - op_slice (sopvar)
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
% DJ, 11/14/2026: Initial coding;

% Slicing needs a 1 x 2 'dims'
dims = P.dims;
if numel(dims)~=2
    error('sopvar:sliceDims',['Slicing needs a 1 x 2 ''dims''; this sopvar '...
        'has %d entries.'],numel(dims))
end
% Linear indices: only a one-column or one-row operator has a sub-operator
% for them. P(:) of a 1 x p operator would be a p x 1 reshape, not a slice.
if isscalar(varargin)
    if dims(2)==1
        varargin = {varargin{1},1};
    elseif dims(1)==1 && ~ischar(varargin{1})
        varargin = {1,varargin{1}};
    else
        error('sopvar:linearIndex',['Linear indexing P(I) needs a sopvar with '...
            'one column, or one row and I not '':''; use P(I,J).'])
    end
elseif numel(varargin)~=2
    error('sopvar:numIndices','Index a sopvar as P(I,J), or P(I) for one row or column.')
end
indr = parse_index(varargin{1},dims(1),'Row');
indc = parse_index(varargin{2},dims(2),'Column');

% Every component in order: P itself, no parameter copied.
if isequal(indr,(1:dims(1))') && isequal(indc,(1:dims(2))')
    Psop = P;
    return
end

% Rows of each parameter for the kept output components and columns for the
% kept input components: component outer, monomial inner.
nL = prod(cellfun(@numel,P.ZL));
nR = prod(cellfun(@numel,P.ZR));
rows = reshape((indr.'-1)*nL + (1:nL).',[],1);
cols = reshape((indc.'-1)*nR + (1:nR).',[],1);

params = P.params;
for ii = 1:numel(params)
    % Full-size cells only; the [] / scalar-0 shorthand stays a zero cell.
    if isequal(size(params{ii}),[dims(1)*nL,dims(2)*nR])
        params{ii} = params{ii}(rows,cols);
    end
end
% Assigned directly: dropping whole components keeps the variable order and
% the canonical multiplier form, so the constructor's checks are not needed.
Psop = P;
Psop.params = params;
Psop.dims = [numel(indr),numel(indc)];

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
if ischar(ind) && isequal(ind,':')
    ind = (1:n)';
    return
end
if islogical(ind),  ind = find(ind);    end
% 'all' of the conditions, so an empty index passes and NaN fails; a
% fractional index would select a misaligned window of rows.
if ~isnumeric(ind) || ~all(ind(:)>=1 & ind(:)<=n & mod(ind(:),1)==0)
    error('sopvar:badIndex',['%s indices into a sopvar must be '':'', '...
        'logical, or integers in 1..%d.'],side,n)
end
ind = double(ind(:));
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
