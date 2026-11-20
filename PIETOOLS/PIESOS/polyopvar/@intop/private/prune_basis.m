function Kop = prune_basis(Kop)
% KOP = PRUNE_BASIS(KOP) removes from the monomials Kop.Z.in those that are
% not used by any of the parameters, in a way that preserves the tensor
% product structure: a monomial of variable i is removed only if none of
% the monomials in the product with the monomials of the other variables
% has a nonzero coefficient. At least one monomial is retained for each
% variable.
%
% INPUTS
% - Kop:    'intop' object;
%
% OUTPUTS
% - Kop:    'intop' object representing the same functional as the input,
%           with the unused monomials removed from Kop.Z.in, and the
%           columns of the parameters Kop.params relabeled accordingly;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - prune_basis
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
% DJ, 10/15/2026: Initial coding;

N = numel(Kop.Z.in);
if N==0
    return
end
[nmons,strides,nZ] = Zdims(Kop.Z.in);
n = size(Kop.params{1},2)/nZ;

% Determine which monomials are used
used = false(nZ,1);
for k=1:numel(Kop.params)
    colnz = full(any(Kop.params{k},1));
    used = used | any(reshape(colnz,nZ,n),2);
end

% Determine which monomials are used along each variable
keep = cell(1,N);
newidx = cell(1,N);
nmons_new = zeros(1,N);
lin = (0:nZ-1)';
zi_all = cell(1,N);
for i=1:N
    zi_all{i} = mod(floor(lin/strides(i)),nmons(i))+1;
    keep{i} = accumarray(zi_all{i},used,[nmons(i),1],@any);
    if ~any(keep{i})
        keep{i}(1) = true;
    end
    newidx{i} = cumsum(keep{i}).*keep{i};
    nmons_new(i) = nnz(keep{i});
end
if all(nmons_new==nmons)
    return
end

% Relabel the columns
nZ_new = prod(nmons_new);
stNew = ones(1,N);
for i=N-1:-1:1
    stNew(i) = stNew(i+1)*nmons_new(i+1);
end
zmap = ones(nZ,1);
for i=1:N
    zmap = zmap + (max(newidx{i}(zi_all{i}),1)-1)*stNew(i);
end
for k=1:numel(Kop.params)
    Kop.params{k} = remap_cols(Kop.params{k},zmap,nZ_new);
end
for i=1:N
    Kop.Z.in{i} = Kop.Z.in{i}(keep{i},:);
end

end
