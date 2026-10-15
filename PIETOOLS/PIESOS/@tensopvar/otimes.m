function C = otimes(A,B,type)
% C = OTIMES(A,B,TYPE) returns the tensor product of the tensor-PI
% operators A and B, so that
%   (A*Zd(x))(s).*(B*Zd(y))(s) = (C*(Zd(x) o Zd(y)))(s).
%
% INPUTS
% - A:      m1 x n1 'tensopvar' object (or 'sopvar', or 1);
% - B:      m2 x n2 'tensopvar' object (or 'sopvar', or 1);
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
% - C:      'tensopvar' object representing the tensor product of A and B;
%
% NOTES
% - If A or B is the scalar 1, the other operator is returned. If A or B
%   has size 0 x 0, C = zeros(0,0) is returned.
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
    type = [1,2];
elseif isscalar(type)
    % Assume only type of row product is specified
    type = [type,2];
end

% Check that the inputs are appropriate
if ~isa(A,'tensopvar')
    if isequal(A,1)
        C = B;
        return
    elseif isa(A,'sopvar')
        A = sopvar2tensopvar(A);
    else
        error("Inputs must be of type 'tensopvar' or 'sopvar'.")
    end
end
if ~isa(B,'tensopvar')
    if isequal(B,1)
        C = A;
        return
    elseif isa(B,'sopvar')
        B = sopvar2tensopvar(B);
    else
        error("Inputs must be of type 'tensopvar' or 'sopvar'.")
    end
end
if ~A.chkval
    A = initialize(A);
end
if ~B.chkval
    B = initialize(B);
end

% Check that the dimensions are appropriate for elementwise multiplication
if type(1)==1 && A.matdim(1)~=B.matdim(1)
    error("Row dimensions of operators must match for elementwise multiplication.")
elseif type(2)==1 && A.matdim(2)~=B.matdim(2)
    error("Column dimensions of operators must match for elementwise multiplication.")
end

% For tensor-PI operators defined by a single factor, there's no
% distinction between pointwise or Kronecker product
if size(A.ops,2)==1
    A.otype = type;
end
if size(B.ops,2)==1
    B.otype = type;
end

% Express the operators in terms of the same set of variables
[A,B] = common_vars(A,B);

C = A;
if all(type==A.otype) && all(type==B.otype)
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
    C.depmat.out = [A.depmat.out; B.depmat.out];
    C.depmat.in = [A.depmat.in; B.depmat.in];
    C.dims = [A.dims;B.dims];
else
    % The products are of a different type
    % --> store the two operators as separate factors
    C.ops = {A,B};
    Adeps = A.deps;
    Bdeps = B.deps;
    if type(1)==1
        C.depmat.out = max(Adeps.out,Bdeps.out);
    else
        C.depmat.out = Adeps.out + Bdeps.out;
    end
    if type(2)==1
        C.depmat.in = max(Adeps.in,Bdeps.in);
    else
        C.depmat.in = Adeps.in + Bdeps.in;
    end
    C.dims = [A.matdim; B.matdim];
    C.otype = type;
    C.chkval = true;
end


end