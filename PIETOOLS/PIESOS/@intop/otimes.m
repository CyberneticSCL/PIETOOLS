function C = otimes(A,B,type)
% C = OTIMES(A,B,TYPE) returns the tensor product of the functionals A
% and B, so that
%   (A*Z_d(x)).*(B*Z_d(y)) = C*(Z_d(x) o Z_d(y))
% where Z_d(x) is a distributed monomial of degree d in the state x (not to
% be confused with the property Zd of 'intop' objects), and o is the
% Kronecker product.
%
% INPUTS
% - A:      m1 x n1 'intop' object, or constant of type 'double',
%           'polynomial', or 'dpvar';
% - B:      m2 x n2 'intop' object, or constant of type 'double',
%           'polynomial', or 'dpvar';
% - type:   (optional) 1 x 2 array of integers indicating what type of
%           product is taken along the row and column dimensions, with 1
%           for an elementwise product and 2 for a Kronecker product. A
%           scalar is read as type(1), with type(2)=2. Defaults to [1,2].
%           If type(1)=1, it is required that m1=m2, or m1=1 or m2=1.
%           Only type(2)=2 is supported;
%
% OUTPUTS
% - C:      m x n1*n2 'intop' object representing the tensor product of A
%           and B, with m=max(m1,m2) if type(1)=1 and m=m1*m2 if
%           type(1)=2;
%
% NOTES
% - If A.vars.in and B.vars.in share a variable, then this is taken to be
%   the same spatial variable in both, with the same domain. The functional
%   C integrates over the dummy variables of both, where the dummy
%   variables of B are numbered after those of A.
% - A constant is treated as a kernel that does not depend on the dummy
%   variables.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - otimes
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

% By default, perform Kronecker product along columns but pointwise product
% along rows
if nargin<2
    error("Not enough inputs")
elseif nargin<=2
    type = [1,2];
elseif isscalar(type)
    % Assume only type of row product is specified
    type = [type,2];
end
if type(2)==1
    error("Elementwise multiplication along the input dimension is not supported.")
end

% Check that both operators are properly declared
if (isa(A,'intop') && ~A.chkval && ~isvalid(A)) || (isa(B,'intop') && ~B.chkval && ~isvalid(B))
    error("One of the operators is not properly declared.")
end

% Make sure the dimensions of the functionals match for elementwise product
m1 = size(A,1);
m2 = size(B,1);
if type(1)==1 && m1~=m2 && m1~=1 && m2~=1
    error("Row dimension of the functionals must match for elementwise product.")
end

% For multiplication with degree-0 monomial, we can just take the product of
% the constant with the kernels, which do not depend on the dummy variables
% of the constant
if isa(A,'double') || isa(A,'polynomial') || isa(A,'dpvar')
    paramsB = kernels(B);
    params = cellfun(@(b) otimes_params(A,b,type),paramsB,'UniformOutput',false);
    C = intop(params,B.intmat,B.vars,B.dom);
    return
elseif ~isa(A,'intop')
    error("Functionals must be specified as 'intop' objects.")
end
if isa(B,'double') || isa(B,'polynomial') || isa(B,'dpvar')
    paramsA = kernels(A);
    params = cellfun(@(a) otimes_params(a,B,type),paramsA,'UniformOutput',false);
    C = intop(params,A.intmat,A.vars,A.dom);
    return
elseif ~isa(B,'intop')
    error("Functionals must be specified as 'intop' objects.")
end

% Establish a unique set of variables in terms of which the product is
% defined
varsA = A.vars.in;          domA = A.dom.in;            NA = numel(varsA); 
varsB = B.vars.in;          domB = B.dom.in;            NB = numel(varsB);
vars_full = [varsA, varsB];   dom_full = [domA;domB];
[vars,idcs_rtn,idcs] = unique(vars_full);
N = numel(vars);
dom = dom_full(idcs_rtn,:);
if any(any(dom(idcs,:)~=dom_full))
    error("Identical variables in functionals must have the same domain.")
end
idcs1 = idcs(1:numel(varsA));       idcs2 = idcs(numel(varsA)+1:end);

% Update parameters to match full variable list
% The kernels are expressed as polynomial/dpvar objects, to take their product
paramsA = kernels(A);               paramsB = kernels(B);
orderA = zeros(1,N+1);   orderA(idcs1) = 1:NA;   orderA(orderA==0) = (NA+1):(N+1);
orderB = zeros(1,N+1);   orderB(idcs2) = 1:NB;   orderB(orderB==0) = (NB+1):(N+1);
paramsA = permute(reshape(paramsA,[size(paramsA,1:NA),ones(1,N+1-NA)]),orderA);
paramsB = permute(reshape(paramsB,[size(paramsB,1:NB),ones(1,N+1-NB)]),orderB);

% Update variables in B to match full variable list
intmatA = repmat({zeros(1,0)},1,N+1);      intmatA(idcs1) = A.intmat;           
intmatB = repmat({zeros(1,0)},1,N+1);      intmatB(idcs2) = B.intmat;
varnames_B_old = cell(0,1);
varnames_B_new = cell(0,1);
% Tensor product in the dummy variables
for i=1:N
    % Update intmats to match full list of variables
    intmatA_i = intmatA{i};     dA = size(intmatA_i,2)/2;
    intmatB_i = intmatB{i};     dB = size(intmatB_i,2)/2;
    intmatA_i(intmatA_i==dA+1) = dA+dB+1;   % index of upper boundary of domain
    intmatB_i = intmatB_i + dA;
    intmatB_i(intmatB_i==dA) = 0;           % index of lower boundary of domain
    intmatA{i} = intmatA_i;
    intmatB{i} = intmatB_i;

    % Determine new variable names for the kernels of B
    for j=1:dB
        varnames_B_old = [varnames_B_old; {[vars{i},'_dum',num2str(j)]}];
        varnames_B_new = [varnames_B_new; {[vars{i},'_dum',num2str(dA+j)]}];
    end
end

% Update dummy variable names in the kernels defining B
for k=1:numel(paramsB)
    if isa(paramsB{k},'double')
        continue
    end
    [~,sub_idcs1,sub_idcs2] = intersect(varnames_B_old,paramsB{k}.varname,'stable');
    new_vars_k = varnames_B_new(sub_idcs1);
    paramsB{k}.varname(sub_idcs2) = new_vars_k;
end

% Combine intmats to consider all combinations of intmatA and intmatB
intmat = cell(1,N);
nterms = ones(1,N);
ntermsA = ones(1,N);
ntermsB = ones(1,N);
for i=1:N
    ntermsA(i) = size(intmatA{i},1);
    ntermsB(i) = size(intmatB{i},1);
    intmat{i} = [repelem(intmatA{i},ntermsB(i),1), repmat(intmatB{i},ntermsA(i),1)];
    nterms(i) = size(intmat{i},1);
end

% Compute kernel for each of the combined integrals
params = cell([nterms,1]);
for k=1:numel(params)
    % Establish which parameter in A and B is associated with this term
    idcs_k = cell(1,N+1);
    [idcs_k{:}] = ind2sub([nterms,1],k);
    idcs_kA = num2cell(ones(1,N+1));
    idcs_kB = num2cell(ones(1,N+1));
    for i=1:N
        [idx_iB,idx_iA] = ind2sub([ntermsB(i),ntermsA(i)],idcs_k{i});
        idcs_kA{i} = idx_iA;
        idcs_kB{i} = idx_iB;
    end
    idxA = sub2ind([ntermsA,1],idcs_kA{:});
    idxB = sub2ind([ntermsB,1],idcs_kB{:});
    % Multiply the specified parameters
    params{k} = otimes_params(paramsA{idxA},paramsB{idxB},type);
end

% Declare the operator representing the tensor product
C = intop(params,intmat,vars,dom,true);

end
