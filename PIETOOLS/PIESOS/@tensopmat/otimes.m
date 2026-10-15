function C = otimes(A,B,type)
% C = OTIMES(A,B,TYPE) returns the tensor product of tensor-PI operators A
% and B, so that (A*Zd(x))(s).*(B*Zd(y))(s) = (C*(Zd(x) o Zd(y)))(s).
%
% INPUTS
% - A:      m1 x n1 'tensopmat' object (or 'tensopvar', 'intop', 'sopvar',
%           or scalar 'double', which multiplies B);
% - B:      m2 x n2 'tensopmat' object (or 'tensopvar', 'intop', 'sopvar',
%           or scalar 'double', which multiplies A);
% - type:   (optional) 1 x 2 array of values 1 or 2 indicating what type of
%           product is taken along the row and column dimensions. Defaults
%           to [1,2]; a scalar specifies type(1), with type(2)=2.
%               If type(1)=1, then a pointwise product is taken along the
%               row dimension, requiring m1=m2. In this case, if A and B
%               share any output variable, then a pointwise product is
%               taken along this variable, i.e.
%                   (Cop*v^2)(s) = (Aop*v)(s)*(Bop*v)(s).
%               If type(1)=2, then a Kronecker product is taken along the
%               row dimension. In this case, if A and B share any output
%               variable, then a tensor product is taken along this
%               spatial direction, i.e.
%                   (Cop*v^2)(s1,s2) = (Aop*v)(s1)*(Bop*v)(s2).
%               Likewise, type(2) specifies the product along the column
%               dimension, for the input variables;
%
% OUTPUTS
% - C:      'tensopmat' object representing the tensor product of A and B;
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

% Check that the inputs are appropriate
if ~isa(A,'tensopmat')
    if isa(A,'double') && isscalar(A)
        % Allow for multiplication with scalar
        C = A*B;
        return
    elseif ~isa(A,'sopvar') && ~isa(A,'tensopvar') && ~isa(A,'intop')
        error("Inputs must be of type 'tensopmat' or 'tensopvar'.")
    else
        A = tensopmat(A);
    end
end
if ~isa(B,'tensopmat')
    if isa(B,'double') && isscalar(B)
        % Allow for multiplication with scalar
        C = A*B;
        return
    elseif ~isa(B,'sopvar') && ~isa(B,'tensopvar') && ~isa(B,'intop')
        error("Inputs must be of type 'tensopmat' or 'tensopvar'.")
    else
        B = tensopmat(B);
    end
end
if ~A.chkval && ~isvalid(A)
    error("Operation not supported: The first 'tensopmat' object is not properly declared.")
else
    A.chkval = true;
end
if ~B.chkval && ~isvalid(B)
    error("Operation not supported: The second 'tensopmat' object is not properly declared.")
else
    B.chkval = true;
end

% Check that the dimensions are appropriate for elementwise multiplication
if type(1)==1 && ~isequal(size(A,1),size(B,1))
    error("Row dimensions of operators must match for elementwise multiplication.")
elseif type(2)==1 && ~isequal(size(A,2),size(B,2)) 
    error("Column dimensions of operators must match for elementwise multiplication.")
end

% Express the operators in terms of the same variables
[A,B] = common_vars(A,B);

% % Compute the appropraite tensor product of each of the operators in the
% % structure
[nrA,ncA] = size(A.ops);
[nrB,ncB] = size(B.ops);
if all(type==2)
    % Proper tensor product along both row and column dimension
    C = A;
    Cops = cell(nrA*nrB,ncA*ncB);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B
        if isempty(A.ops{idx})
            continue
        end
        [ridx,cidx] = ind2sub(size(A.ops),idx);
        ABops = cellfun(@(a,rb) otimes_op(A.ops{idx},a,type,A,ridx,B,rb),B.ops,num2cell(repmat((1:nrB)',1,ncB)),'UniformOutput',false);
        Cops((ridx-1)*nrB+(1:nrB),(cidx-1)*ncB+(1:ncB)) = ABops;
    end
    C.ops = Cops;
    % Set the dimensions
    C.dims.in = repelem(A.dims.in,ncB,1).*repmat(B.dims.in,ncA,1);
    C.dims.out = repelem(A.dims.out,nrB,1).*repmat(B.dims.out,nrA,1);
    % Set the dependency arrays
    C.depmat.in = repelem(A.depmat.in,ncB,1) + repmat(B.depmat.in,ncA,1);
    C.depmat.out = repelem(A.depmat.out,nrB,1) + repmat(B.depmat.out,nrA,1);
    C.chkval = true;

elseif type(1)==2
    % Pointwise product along column dimension
    C = A;
    Cops = cell(nrA*nrB,ncA);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B(:,j)
        if isempty(A.ops{idx})
            continue
        end
        [ridx,cidx] = ind2sub(size(A.ops),idx);
        ABops = cellfun(@(a,rb) otimes_op(A.ops{idx},a,type,A,ridx,B,rb),B.ops(:,cidx),num2cell((1:nrB)'),'UniformOutput',false);
        Cops((ridx-1)*nrB+(1:nrB),cidx) = ABops;
    end
    C.ops = Cops;
    % Set the dimensions
    C.dims.in = A.dims.in;
    C.dims.out = repelem(A.dims.out,nrB,1).*repmat(B.dims.out,nrA,1);
    % Set the dependency arrays
    C.depmat.in = max(A.depmat.in,B.depmat.in);
    C.depmat.out = repelem(A.depmat.out,nrB,1) + repmat(B.depmat.out,nrA,1);
    C.chkval = true;
    
elseif type(2)==2
    % Pointwise product along row dimension
    C = A;
    Cops = cell(nrA,ncA*ncB);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B(i,:)
        if isempty(A.ops{idx})
            continue
        end
        [ridx,cidx] = ind2sub(size(A.ops),idx);
        ABops = cellfun(@(a) otimes_op(A.ops{idx},a,type,A,ridx,B,ridx),B.ops(ridx,:),'UniformOutput',false);
        Cops(ridx,(cidx-1)*ncB+(1:ncB)) = ABops;
    end
    C.ops = Cops;
    % Set the dimensions
    C.dims.in = repelem(A.dims.in,ncB,1).*repmat(B.dims.in,ncA,1);
    C.dims.out = A.dims.out;
    % Set the dependency arrays
    C.depmat.in = repelem(A.depmat.in,ncB,1) + repmat(B.depmat.in,ncA,1);
    C.depmat.out = max(A.depmat.out,B.depmat.out);
    C.chkval = true;

else
    % Pointwise product along both row and column dimension
    C = A;
    Cops = cell(nrA,ncA);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B{i,j}
        if isempty(A.ops{idx}) || isempty(B.ops{idx})
            continue
        end
        [ridx,~] = ind2sub(size(A.ops),idx);
        Cops{idx} = otimes_op(A.ops{idx},B.ops{idx},type,A,ridx,B,ridx);
    end
    C.ops = Cops;
    % Set the dimensions
    C.dims.in = A.dims.in;
    C.dims.out = A.dims.out;
    % Set the dependency arrays
    C.depmat.in = max(A.depmat.in,B.depmat.in);
    C.depmat.out = max(A.depmat.out,B.depmat.out);
    C.chkval = true;
end

end

function C = otimes_op(A,B,type,XA,iA,XB,iB)
% Tensor product of two operators, where a scalar constant (a 'tensopvar'
% consisting of a 'sopvar' without variables, representing a degree-0
% monomial) multiplies the kernels of the other operator, as a double would
% (intop/otimes does so too, but checks A.chkval before its double branch)
% A is block iA of the 'tensopmat' XA, B block iB of XB (row indices), for the output space of a
% block of type double/polynomial/dpvar, which is converted to an operator (native2op)
if is_native_op(A) && is_native_op(B)
    C = otimes_native(A,B,type);        % two matrices acting on R
    return
elseif is_native_op(A)
    A = native2op(A,XA,iA);             % multiplier operator, a tensopvar
    if isa(B,'intop') && ~is_constant(A)
        error("The tensor product of a block of type 'double' or 'polynomial' with a functional ('intop') is only supported for a scalar constant.")
    end
elseif is_native_op(B)
    B = native2op(B,XB,iB);
    if isa(A,'intop') && ~is_constant(B)
        error("The tensor product of a block of type 'double' or 'polynomial' with a functional ('intop') is only supported for a scalar constant.")
    end
end
if isa(A,'intop') && is_constant(B)
    c = B.ops{1}.params{1};
    C = A.*c;
elseif isa(B,'intop') && is_constant(A)
    c = A.ops{1}.params{1};
    C = c.*B;
elseif type(1)==1 && is_constant(A) && is_constant(B)
    % The product of two degree-0 monomials is a degree-0 monomial (otimes of
    % two tensopvar would give a product of two factors, applied to a state)
    % Only for the pointwise product along rows, as in the product of polynomials: for
    % a Kronecker product (e.g. composing the factors of a monomial in subs), the two
    % factors must be kept, as one may be a state without spatial variables.
    C = mat2tensopvar(A.ops{1}.params{1}.*B.ops{1}.params{1});
else
    C = otimes(A,B,type);
end
end

function is_const = is_constant(A)
% Check whether A is a scalar constant, as declared by `constant2polyopvar`
is_const = isa(A,'tensopvar') && isscalar(A.ops) && isa(A.ops{1},'sopvar') && ...
            all(size(A)==1) && isempty(A.ops{1}.vars.in) && isempty(A.ops{1}.vars.out) && ...
            isscalar(A.ops{1}.params);
end

function C = otimes_native(A,B,type)
% C = OTIMES_NATIVE(A,B,TYPE) returns the tensor product of two matrices A, B of
% type 'double', 'polynomial' or 'dpvar', acting on R: for type(1)=1 the rows are
% multiplied pointwise (as the product of the functions, requiring equal numbers of
% rows), for type(1)=2 the Kronecker product is taken; the same for the columns
% with type(2). The result is the matrix of products A(r1,c1)*B(r2,c2) in the
% row-column order of the Kronecker product, with r1=r2 (c1=c2) for a pointwise
% product.
%

if type(1)==2 && ~isempty(intersect(native_vars(A),native_vars(B)))
    error("The tensor product along the rows of blocks of type 'double', 'polynomial' or 'dpvar' that share a spatial variable is not supported: use a pointwise product.")
end
[mA,nA] = size(A);      [mB,nB] = size(B);
if type(1)==1
    if mA~=mB
        error("Row dimensions of operators must match for elementwise multiplication.")
    end
    ridcsA = 1:mA;      ridcsB = 1:mB;
else
    ridcsA = repelem(1:mA,mB);      ridcsB = repmat(1:mB,1,mA);
end
if type(2)==1
    if nA~=nB
        error("Column dimensions of operators must match for elementwise multiplication.")
    end
    cidcsA = 1:nA;      cidcsB = 1:nB;
else
    cidcsA = repelem(1:nA,nB);      cidcsB = repmat(1:nB,1,nA);
end
C = A(ridcsA,cidcsA).*B(ridcsB,cidcsB);

end
