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
