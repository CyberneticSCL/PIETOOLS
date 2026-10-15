function V = innerprod(u1,P,u2)
% V = INNERPROD(U1,P,U2) returns the weighted inner product of two
% distributed polynomial functions, V = <U1,P*U2>.
%
% INPUTS
% - u1,u2:  n x 1 'polyopvar' objects representing distributed polynomial
%           functions. If u2 = [] or omitted, u2 = u1 is used;
% - P:      n x n object of type 'double', 'polynomial', or 'dpvar', or a
%           cell array of such blocks, partitioned as the vectors u1 and
%           u2 (block (i,j) has size u1.dims.out(i) x u2.dims.out(j)).
%           A 'polynomial' or 'dpvar' may depend on the spatial
%           variables, and is then integrated along with the functions.
%           A scalar P is used as P times the identity. If only u1 is
%           passed, P = 1. If INNERPROD(U1,U2) is called with a
%           'polyopvar' U2 in the second position, P = 1 and u2 = U2;
%
% OUTPUTS
% - V:      1 x 1 'polyopvar' object representing the inner product of the
%           distributed polynomial vectors with weight P.
%
% NOTES
% - If u1 is not a 'polyopvar', V is computed as (u1'*P)*u2. If u2 is not
%   a 'polyopvar', V is computed as (P*u2)'*u1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - innerprod
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
% DJ, 10/15/2026: Initial coding;

% Check that the inputs are appropriate
is_symmetric = false;
if nargin==1
    u2 = u1;
    P = 1;
    is_symmetric = true;
end
% u2 omitted or empty: take u2 = u1, innerprod(u1,P) = innerprod(u1,P,u1)
if (nargin==2 && ~isa(P,'polyopvar')) || (nargin>=3 && isempty(u2))
    u2 = u1;
    is_symmetric = isnumeric(P);
end
if nargin<=2 && isa(P,'polyopvar')
    % Assume P represents a vector rather than the weight
    u2 = P;
    P = 1;
end
if ~isa(u1,'polyopvar') 
    V = (u1'*P)*u2;
    return
elseif ~isa(u2,'polyopvar')
    V = (P*u2)'*u1;
    return
end
if ~is_symmetric
    [u1,u2] = common_vars(u1,u2);
end
if ~isnumeric(P) && ~isa(P,'polynomial') && ~isa(P,'dpvar') && ~isa(P,'cell')
    error("Inner product weight must be specified as object of type 'double', 'polynomial', or 'dpvar'.")
end
if is_symmetric && ~all(all(isequal(P-P',zeros(size(P)))))
    is_symmetric = false;
end


% Make sure the dimensions are appropriate for taking the inner product
n_blks1 = numel(u1.dims.out);
n_blks2 = numel(u2.dims.out);
if isscalar(P) && ~isa(P,'cell')
    if ~isequal(u1.dims.out,u2.dims.out)
        error("Dimensions of vectors for which to compute inner product must match.")
    elseif ~isequal(u1.vars.out,u2.vars.out) || ~isequal(u1.dom.out,u2.dom.out) || ~isequal(u1.depmat.out,u2.depmat.out)
        error("Output function spaces of the vectors must match.")
    end
else
    if isa(P,"cell")
        if size(P,1)~=n_blks1 || size(P,2)~=n_blks2
            error("Partitioning of the weight must match partioning of the vectors.")
        end
        nr_arr = cellfun(@(a) size(a,1),P);
        nc_arr = cellfun(@(a) size(a,2),P);
        % Block P{i,j} is dims1(i) x dims2(j): rows constant along block 
        % rows, columns along block columns
        if any(any(nr_arr - nr_arr(:,1))) || ~isequal(u1.dims.out,nr_arr(:,1)) 
            error("Row dimensions of the weights in each column of blocks should match dimensions of left vector.")
        elseif any(any(nc_arr - nc_arr(1,:))) || ~isequal(u2.dims.out,nc_arr(1,:)')
            error("Column dimensions of the weights in each row of blocks should match dimensions of right vector.")
        end
    elseif size(P,1)~=size(u1,1) || size(P,2)~=size(u2,1)
        error("Dimensions of the vectors should match dimensions of the weight.")
    end
end


% Split the polynomials into coefficients and a monomial basis
C1 = u1.C;      Z1 = u1.Z;
C2 = u2.C;      Z2 = u2.Z;

% Declare the full list of monomials appearing in the product of Z1 and Z2
[ridcs1,cidcs1,Z1vals] = find(Z1);     [nZ1,nc1] = size(Z1);
[ridcs2,cidcs2,Z2vals] = find(Z2);     [nZ2,nc2] = size(Z2);
% find returns row vectors for a row-vector Z; column form assumed below
ridcs1 = ridcs1(:);  cidcs1 = cidcs1(:);  Z1vals = Z1vals(:);
ridcs2 = ridcs2(:);  cidcs2 = cidcs2(:);  Z2vals = Z2vals(:);
% The factors of u2 follow all factors of u1, including multipliers from R
% (zeros in Z1) 
nfac1 = num_factors(u1); 
d1 = repelem(nfac1,nZ2,1); 
ridcs1 = (ridcs1-1)*nZ2 + (1:nZ2);                  % repelem(Z1,nz2,1)
cidcs1 = repmat(cidcs1,1,nZ2);
Z1vals = repmat(Z1vals,1,nZ2);
ridcs2 = ridcs2 + (0:nZ1-1)*nZ2;                    % repmat(Z2,nz1,1);
% d1(ridcs2) has the orientation of d1 if ridcs2 is a vector (nnz(Z2)==1)
cidcs2 = repmat(cidcs2,1,nZ1) + reshape(d1(ridcs2),size(ridcs2));
Z2vals = repmat(Z2vals,1,nZ1);
Z_full = sparse([ridcs1(:);ridcs2(:)],[cidcs1(:);cidcs2(:)],[Z1vals(:);Z2vals(:)],nZ1*nZ2,max([nc1;nfac1(:)+nc2]));

% Retain only unique monomials
[M,Z_lin] = uniquerows_integerTable([sum(Z_full>0,2),Z_full]);  % Z_full = M*Z_lin
Z_lin = Z_lin(:,2:end);
nZ_lin = size(Z_lin,1);
old2new_idcs = M*(1:nZ_lin)';               % Z_full = Z_lin(old2new_idcs,:);


% Get the dimensions of the different operators
blkdims1 = u1.dims.out;
blkdims1_cum = [0;cumsum(blkdims1)];
blkdims2 = u2.dims.out;
blkdims2_cum = [0;cumsum(blkdims2)];
blk_nums = 1:n_blks1*n_blks2;
if is_symmetric
    blk_nums = reshape(tril(reshape(blk_nums,n_blks1,n_blks2)),1,[]);
    % tril leaves zeros for the upper triangle
    blk_nums(blk_nums==0) = [];
end

% Split the weight into blocks
if isa(P,'cell')
    Pcell = P;
else
    Pcell = cell(n_blks1,n_blks2);
    if isscalar(P)
        % Diagonal blocks P, off-diagonal blocks the double 0 
        % (skipped in the loop below)
        Pcell = num2cell(zeros(n_blks1));
        for i=1:n_blks1
            Pcell{i,i} = P;  
        end
    else
        for k=1:numel(Pcell)
            [i,j] = ind2sub(size(Pcell),k);
            Pcell{i,j} = P(blkdims1_cum(i)+1:blkdims1_cum(i+1), blkdims2_cum(j)+1:blkdims2_cum(j+1));
        end
    end
end

% Initialize an empty coefficient operator acting on the new basis
Kcell = num2cell(zeros(1,nZ_lin));
for k=blk_nums
    % Extract the row and column index of the parameter block
    % Pcell{k} is column-major in (n_blks1,n_blks2)
    [ridx,cidx] = ind2sub([n_blks1,n_blks2],k);

    % Extract the coefficient operators associated with this monomial
    C1_i = C1.ops(ridx,:);
    C2_j = C2.ops(cidx,:);

    % Extract the kth block of decision variables
    Pij = Pcell{k};
    if all(all(isequal(Pij,zeros(size(Pij)))))
        continue
    end

    for l=1:nZ1*nZ2
        % Determine which monomials in u1 and u2 correspond with monomial l
        [z2,z1] = ind2sub([nZ2,nZ1],l);

        % Extract the coefficient operators associated with this monomial
        C1_il = C1_i{z1};
        C2_jl = C2_j{z2};
        if isempty(C1_il) || isempty(C2_jl)
            continue
        end

        % Take the inner product of the coefficient operators
        if isscalar(Pij)
            % Weight is passed to int, so that it is integrated if it 
            % varies in space
            Kij = int(otimes(C1_il,C2_jl,[1,2]),Pij*ones(1,blkdims1(ridx)));
        else
            m = size(C1_il,1);   n = size(C2_jl,1);
            perm = reshape(reshape(1:(m*n),m,n).',1,[]); 
            % Pij(perm) is a column if Pij has one column; Kij needs a row
            Pperm = Pij(perm);      
            if size(Pperm,1)>1
                Pperm = Pperm.';
            end
            Kij = int(otimes(repelem(C1_il,n,1),repmat(C2_jl,m,1),[1,2]),Pperm);
        end

        % Account for symmetry
        if is_symmetric && ridx~=cidx
            Kij = 2*Kij;
        end

        % Match the obtained kernels with the correct distributed monomial
        mon_idx = old2new_idcs(l);
        Kcell{mon_idx} = Kcell{mon_idx} + Kij;
    end    
end

% Declare the coefficient matrix defining V
% Drop monomials that received no contribution (still the double 0), e.g. 
% from zero off-diagonal blocks of the weight, which tensopmat cannot 
% convert
is_nonzero = ~cellfun(@isnumeric,Kcell);
Kcell = Kcell(is_nonzero);      Z_lin = Z_lin(is_nonzero,:);
% No monomial contributes, i.e. V is identically zero
if isempty(Kcell)
    V = constant2polyopvar(0);
    return
end
Kop = tensopmat(Kcell);

% Declare the polynomial defined by coefficients K and monomials Z
V = u1;
V.Z = Z_lin;
V.C = Kop;
V.dims.out = 1;
V.vars.out = Kop.vars.out;
V.dom.out = Kop.dom.out;
V.depmat.out = Kop.depmat.out;
V.chkval = true;

% Get rid of duplicate terms in the functional
V = combine_terms(V);

end