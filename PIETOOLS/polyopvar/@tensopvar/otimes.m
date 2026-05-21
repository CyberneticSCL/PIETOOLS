function C = otimes(A,B,type)
% C = OTIMES(A,B) returns the tensor product of coefficient operators A and
% B, so that (A*Zd(x))(s).*(B*Zd(y))(s) = (C*(Zd(x) o Zd(y)))(s)
%
% INPUTS
% - A:  m1 x n1 'tensopvar' object
% - B:  m2 x n2 'tensopvar' object
% - type:   1 x 2 logical array indicating what type of product is taken
%           along the row and column dimensions. 
%               If type(1)=0, then a pointwise product is taken along the 
%               row dimension, requiring m1=m2. In this case, if 
%               A.vars(:,1) and B.vars(:,2) share any variable, then a 
%               pointwise product is taken along this variable, i.e. 
%                   (Cop*v^2)(s) = (Aop*v)(s)*(Bop*v)(s).
%               If type(1)=1, then a Kronecker product is taken along the
%               row dimension. In this case, if A.vars(:,1) and B.vars(:,2) 
%               share any variable, then a tensor product is taken along
%               this spatial direction, i.e.
%                   (Cop*v^2)(s1,s2) = (Aop*v)(s1)*(Bop*v)(s2).
%
% OUTPUTS
% - C:  'tensopvar' object representing the tensor product of A and B. 
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
% DJ, 04/15/2026: Initial coding

% Allow for tensor product with empty array
if nargin<2
    error("Not enough inputs")
elseif all(size(A)==0) || all(size(B)==0)
    C = zeros(0,0);
    return
end

% By default, perform Kronecker product along columns but pointwise product
% along rows
if nargin<=2
    type = [false,true];
elseif isscalar(type)
    % Assume only type of row product is specified
    type = [type,true];
end

% Check that the inputs are appropriate
if ~isa(A,'tensopvar')
    if isequal(A,1)
        C = B;
        return
    elseif ~isa(A,'nopvar') && ~isa(A,'ndopvar')
        error("Inputs must be of type 'tensopvar' or 'nopvar'.")
    else
        A = ndopvar2tensopvar(A);
    end
elseif ~isa(B,'tensopvar')
    if isequal(B,1)
        C = A;
        return
    elseif ~isa(B,'nopvar') && ~isa(B,'ndopvar')
        error("Inputs must be of type 'tensopvar' or 'nopvar'.")
    else
        B = ndopvar2tensopvar(B);
    end
end

% Check that the dimensions are appropriate for elementwise multiplication
if ~type(1) && A.dim(1)~=B.dim(1)
    error("Row dimensions of operators must match for elementwise multiplication.")
elseif ~type(2) && A.dim(2)~=B.dim(2)
    error("Column dimensions of operators must match for elementwise multiplication.")
end

% Express the two operators in terms of the same spatial and dummy
% variables
[A,B] = common_vars(A,B);

% For tensor-PI operators defined by a single factor, there's no
% distinction between pointwise or Kronecker product
if size(A.ops,2)==1
    A.type = type;
end
if size(B.ops,2)==1
    B.type = type;
end

C = A;
if all(type==A.type) && all(type==B.type)
    % The proposed tensor product is of the same type in the definition of
    % A and B
    % --> we just concatenate the factors to represent their product
    [nr_A,nc_A] = size(A.ops);
    [nr_B,nc_B] = size(B.ops);
    C.ops = cell(nr_A*nr_B,nc_A+nc_B);
    for i=1:nr_A*nr_B
        [i_A,i_B] = ind2sub([nr_A,nr_B],i);
        C.ops(i,:) = [A.ops(i_A,:),B.ops(i_B,:)];
    end
    % Keep track of which variables each factor depends on
    C.depmat1 = [A.depmat1; B.depmat1];
    C.depmat2 = [A.depmat2; B.depmat2];
    % Set the order in which the factors should appear
    C.order = [A.order, B.order+nc_A];
else
    % The products are of a different type
    % --> store the two operators as separate factors
    C.ops = {A,B};
    A_depmat = A.dep;
    B_depmat = B.dep;
    C.depmat1 = [A_depmat(1,:),B_depmat(1,:)];
    C.depmat2 = [A_depmat(2,:),B_depmat(2,:)];
    C.type = type;
    C.order = [1,2];
end


end