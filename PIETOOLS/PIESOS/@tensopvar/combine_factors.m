function C = combine_factors(C)
% C = COMBINE_FACTORS(C) merges the factors of the 'tensopvar' object C
% that are operators from R^n to R^m (without spatial variables) into a
% neighboring factor, e.g. for the tensor product A o B, with A an
% operator from R^n to R^m,
%   (A o B)(x o y) = kron(A,B)(x o y),
% with kron(A,B) the operator whose kernels are the Kronecker products of A
% with the kernels of B. This reduces the number of factors of C.
%
% INPUTS
% - C:      'tensopvar' object representing a tensor-PI operator
%               C = sum_{i=1}^{K} Ri1 \otimes ... \otimes Rid,
%           with a proper tensor product of the inputs (C.otype(2)=2), and
%           either a Kronecker product (C.otype(1)=2) or a pointwise
%           product (C.otype(1)=1) of the outputs;
%
% OUTPUTS
% - C:      'tensopvar' object representing the same operator, but with
%           each factor that maps from R^n to R^m, with no input or output
%           variables, merged with a neighboring 'sopvar' factor, the
%           factor to its right, or to its left if it is the last factor.
%           The merged factor acts on the Kronecker product of the inputs
%           of the original factors. Its output is the Kronecker product
%           of their outputs if C.otype(1)=2, and the elementwise product
%           if C.otype(1)=1, for which the factors have the same number of
%           rows, so that row r of the merged operator is the Kronecker
%           product of the rows r of the two operators. A factor from R is
%           merged in every term of C;
%
% NOTES
% - If all factors map from and to R, they are merged into a single factor.
% - C is returned as it is if C.otype(2)=1, if C has a single factor, or if
%   no factor maps from and to R. An elementwise product of the inputs,
%   (A z).*(B z) for the same input z, is quadratic in z, and cannot be
%   expressed as a single linear operator. A factor that maps from and to R
%   is also left as it is if its neighbors are not 'sopvar' objects with
%   numeric parameters, e.g. nested 'tensopvar' factors or 'intop' objects
%   with decision variables, and so are those neighbors.
% - The merged factor has the variables, domains and monomial bases of the
%   neighbor, so that the cost is that of the Kronecker product of the
%   (sparse) parameters, linear in their number of nonzeros: it does not
%   depend on the number of decision variables, which 'sopvar' objects do
%   not have.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - combine_factors (tensopvar)
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
% DJ, 11/08/2026: Initial coding;

if ~isa(C,'tensopvar')
    error("Input must be of type 'tensopvar'.")
end
if ~C.chkval
    C = initialize(C);
end

% The inputs must be combined by a proper tensor product, for which the
% mixed product property gives the Kronecker structure of the merged
% operator. For a pointwise product of the inputs, the operator is not
% linear in the input
if C.otype(2)~=2 || size(C.ops,2)<=1
    return
end

% Determine which factors are 'sopvar' objects with numeric parameters in
% every term, and which of these map from and to R (no input or output
% variables)
is_sop = false(1,size(C.ops,2));
is_R = false(1,size(C.ops,2));
for k=1:size(C.ops,2)
    is_sop(k) = all(cellfun(@(c) isa(c,'sopvar') && ...
                    all(cellfun(@(p) isnumeric(p) || islogical(p),c.params(:))),C.ops(:,k)));
    is_R(k) = is_sop(k) && all(cellfun(@(c) isempty(c.vars.in) && isempty(c.vars.out),C.ops(:,k)));
end

% Merge each factor from R into its neighbor, the next factor if it exists
is_changed = false;
k = 1;
while k<=size(C.ops,2) && size(C.ops,2)>1
    if ~is_R(k)
        k = k+1;
        continue
    end
    if k<size(C.ops,2) && is_sop(k+1)
        j = k+1;
    elseif k>1 && is_sop(k-1)
        j = k-1;
    else
        % No neighbor that can be merged with
        k = k+1;
        continue
    end
    for i=1:size(C.ops,1)
        C.ops{i,j} = merge_sopvar(C.ops{i,k},C.ops{i,j},j>k,C.otype(1));
    end
    % The merged factor maps from R if both did
    is_R(j) = is_R(j) && is_R(k);
    C.ops(:,k) = [];        is_R(k) = [];       is_sop(k) = [];
    is_changed = true;
    % Check the merged factor again, if it is still from R
    if j<k
        k = k-1;
    end
end

% Update the dimensions, variables and dependence of the factors
if is_changed
    C = initialize(C);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function P = merge_sopvar(R,B,R_first,otype)
% P = MERGE_SOPVAR(R,B,R_FIRST,OTYPE) returns the 'sopvar' object P
% representing the product of the operator R from R^n to R^m with the
% operator B, so that P(x o y) = (R x) o (B y) if R_FIRST, and
% P(y o x) = (B y) o (R x) otherwise, for x in R^n and y in the domain of
% B, where o is the Kronecker product of the inputs, and of the outputs
% unless OTYPE==1, in which case the outputs are multiplied elementwise.
%
% INPUT
% - R:      'sopvar' object with no input or output variables, with a
%           single parameter, a numeric m x n matrix;
% - B:      'sopvar' object with numeric parameters, with m rows if
%           otype==1;
% - R_first: logical, true if R is the first factor of the product, i.e.
%           the outer one in the Kronecker product;
% - otype:  integer 1 for a pointwise product and 2 for a tensor product
%           of the outputs;
%
% OUTPUT
% - P:      'sopvar' object with the variables, domains and monomials of B,
%           and with dimensions [m*B.dims(1), n*B.dims(2)], or
%           [B.dims(1), n*B.dims(2)] if OTYPE==1. The parameters are
%           those of B, with each entry multiplied by an entry of R.
%           Rows and columns of the parameters are ordered with the
%           component outer and the monomial inner, so that for R_FIRST
%           the parameters are kron(R,Pg), or, if OTYPE==1, the row-wise
%           Kronecker products of the rows of R and Pg. Otherwise, the
%           index of R is moved between the components of B and the
%           monomials, by remapping the nonzero entries;

mR = R.dims(1);     nR = R.dims(2);
Rmat = R.params{1};
if isempty(Rmat) || (isscalar(Rmat) && Rmat==0)
    Rmat = sparse(mR,nR);                 % zero parameter
else
    Rmat = sparse(Rmat);
end
is_pointwise = otype(1)==1;
if is_pointwise && mR~=B.dims(1)
    error("Row dimensions of the factors must match for an elementwise product.")
end

% Numbers of monomials in the rows and columns of the parameters of B
nL = prod(cellfun(@numel,B.ZL));
nR_B = prod(cellfun(@numel,B.ZR));

params = B.params;
for g=1:numel(params)
    Pg = params{g};
    if isempty(Pg) || (isscalar(Pg) && Pg==0)
        continue                        % zero parameter
    end
    if is_pointwise
        % Entry (r,c) of Pg, in row r of the component ceil(r/nL), is
        % multiplied with the entries (ceil(r/nL),j) of R, moving to the
        % column ((j,k),mono) if R comes first, and to ((k,j),mono)
        % otherwise, for the column (k,mono) of Pg
        [ridcs,cidcs,vals] = find(sparse(Pg));
        [t,j,vR] = find(Rmat(ceil(ridcs/nL),:));
        if R_first
            cidcs_new = (j-1)*size(Pg,2) + cidcs(t);
        else
            k = floor((cidcs(t)-1)/nR_B);       mm = mod(cidcs(t)-1,nR_B);
            cidcs_new = (k*nR + (j-1))*nR_B + mm + 1;
        end
        params{g} = sparse(ridcs(t),cidcs_new,vR.*vals(t),size(Pg,1),nR*size(Pg,2));
    elseif R_first
        params{g} = kron(Rmat,sparse(Pg));
    else
        % Rows ((iB,m),iR) of kron(Pg,R), for the component iB of B,
        % monomial m and component iR of R, move to ((iB,iR),m), and the
        % same for the columns
        [ridcs,cidcs,vals] = find(kron(sparse(Pg),Rmat));
        t = floor((ridcs-1)/mR);        iR = mod(ridcs-1,mR);
        iB = floor(t/nL);               m = mod(t,nL);
        r_new = (iB*mR + iR)*nL + m + 1;
        u = floor((cidcs-1)/nR);        jR = mod(cidcs-1,nR);
        jB = floor(u/nR_B);             mm = mod(u,nR_B);
        cidcs_new = (jB*nR + jR)*nR_B + mm + 1;
        params{g} = sparse(r_new,cidcs_new,vals,mR*B.dims(1)*nL,nR*B.dims(2)*nR_B);
    end
end

if is_pointwise
    P = sopvar(params,B.vars,B.ZL,B.ZR,B.dom,[B.dims(1),nR*B.dims(2)]);
else
    P = sopvar(params,B.vars,B.ZL,B.ZR,B.dom,[mR*B.dims(1),nR*B.dims(2)]);
end

end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
