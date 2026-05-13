function C = otimes(A,B,type)
% C = OTIMES(A,B) returns the tensor product of tensor-PI operators A and
% B, so that (A*Zd(x))(s).*(B*Zd(y))(s) = (C*(Zd(x) o Zd(y)))(s)
%
% INPUTS
% - A:  m1 x n1 'tensopmat' object
% - B:  m2 x n2 'tensopmat' object
% - type:   1 x 2 logical array indicating what type of product is taken
%           along the row and column dimensions. 
%               If type(1)=0, then a pointwise product is taken along the 
%               row dimension, requiring m1=m2. In this case, if 
%               A.vars(:,1) and B.vars(:,2) share any variable, then a 
%               pointwise product is taken along this variable, i.e. 
%                   (Cop*v^2)(s) = (Aop*v)(s)*(Bop*v)(s).
%               If type(1)=1, then a kronecker product is taken along the
%               row dimension. In this case, if A.vars(:,1) and B.vars(:,2) 
%               share any variable, then a tensor product is taken along
%               this spatial direction, i.e.
%                   (Cop*v^2)(s1,s2) = (Aop*v)(s1)*(Bop*v)(s2).
%
% OUTPUTS
% - C:  'tensopmat' object representing the tensor product of A and B. 
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
% DJ, 04/22/2026: Initial coding


% By default, perform Kronecker product along columns but pointwise product
% along rows
if nargin<2
    error("Not enough inputs")
elseif nargin<=2
    type = [false,true];
elseif isscalar(type)
    % Assume only type of row product is specified
    type = [type,true];
end

% Check that the inputs are appropriate
if ~isa(A,'tensopmat')
    if ~isa(A,'tensopvar') && ~isa(A,'intop')
        error("Inputs must be of type 'tensopmat' or 'tensopvar'.")
    else
        A = tensopmat(A);
    end
elseif ~isa(B,'tensopmat')
    if ~isa(B,'tensopvar') && ~isa(B,'intop')
        error("Inputs must be of type 'tensopmat' or 'tensopvar'.")
    else
        B = tensopmat(B);
    end
end

% Check that the dimensions are appropriate for elementwise multiplication
if ~type(1) && ~isequal(A.dim{1},B.dim{1})
    error("Row dimensions of operators must match for elementwise multiplication.")
elseif ~type(2) && ~isequal(A.dim{2},B.dim{2})
    error("Column dimensions of operators must match for elementwise multiplication.")
end

% % Compute the appropraite tensor product of each of the operators in the
% % structure
[nrA,ncA] = size(A.ops);
[nrB,ncB] = size(B.ops);
if type(1) && type(2)
    % Proper tensor product along both row and column dimension
    Cops = cell(nrA*nrB,ncA*ncB);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B
        if isempty(A.ops{idx})
            continue
        end
        [ridx,cidx] = ind2sub(size(A.ops),idx);
        ABops = cellfun(@(a) otimes(A.ops{idx},a,type),B.ops,'UniformOutput',false);
        Cops((ridx-1)*nrB+(1:nrB),(cidx-1)*ncB+(1:ncB)) = ABops;
    end
elseif type(1)
    % Pointwise product along column dimension
    Cops = cell(nrA*nrB,ncA);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B(:,j)
        if isempty(A.ops{idx})
            continue
        end
        [ridx,cidx] = ind2sub(size(A.ops),idx);
        ABops = cellfun(@(a) otimes(A.ops{idx},a,type),B.ops(:,cidx),'UniformOutput',false);
        Cops((ridx-1)*nrB+(1:nrB),cidx) = ABops;
    end
    
elseif type(2)
    % Pointwise product along row dimension
    Cops = cell(nrA,ncA*ncB);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B(i,:)
        if isempty(A.ops{idx})
            continue
        end
        [ridx,cidx] = ind2sub(size(A.ops),idx);
        ABops = cellfun(@(a) otimes(A.ops{idx},a,type),B.ops(ridx,:),'UniformOutput',false);
        Cops(ridx,(cidx-1)*ncB+(1:ncB)) = ABops;
    end
else
    % Pointwise product along both row and column dimension
    Cops = cell(nrA,ncA);
    for idx=1:numel(A.ops)
        % Take tensor product of A{i,j} with B{i,j}
        if isempty(A.ops{idx}) || isempty(B.ops{idx})
            continue
        end
        Cops{idx} = otimes(A.ops{idx},B.ops{idx},type);
    end
end

% Store the tensor-PI operators in a tensopmat object
C = tensopmat();
C.ops = Cops;

% Determine the variables in terms of which the operator is defined
varname1_A = pvar2varname(A.vars(:,1));
varname1_B = pvar2varname(B.vars(:,1));
[~,idcs] = unique([varname1_A; varname1_B],'stable');
vars_C = [A.vars; B.vars];
C.vars = vars_C(idcs,:);
dom_C = [A.dom; B.dom];
C.dom = dom_C(idcs,:);

% Determine the dependency array
[depmat1,depmat2] = get_deps(C);
C.depmat1 = depmat1;
C.depmat2 = depmat2;

end