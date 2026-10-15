function C = times(A,B)
% C = TIMES(A,B) returns the 'intop' object C representing the elementwise
% product of the multiplier operator defined by A with the functional
% defined by B.
%
% INPUTS
% - A:      m x n 'double', 'polynomial', or 'dpvar' object representing
%           some multiplier operator to multiply with B;
% - B:      m x n 'intop' object;
%
% OUTPUTS
% - C:      m x n 'intop' object representing the elementwise product of
%           the multiplier A with the functional defined by B;
%
% NOTES
% - For a constant ('double') multiplier A, the nonzeros of the sparse
%   parameters are multiplied by the entry of A for their matrix row and
%   column. For a 'polynomial' or 'dpvar' multiplier, the kernels are
%   multiplied and expressed in the same monomials and decision variables
%   again.
% - A dimension of size 1 of A or B is repeated in the elementwise product,
%   so that A may also be a scalar. The arguments may be given in either
%   order.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - times
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

% Check that multiplication is supported
if isa(A,'intop') && isa(B,'intop')
    error("Multiplication of functionals is not supported.")
elseif isa(A,'intop')
    % Make sure B is the functional
    C = B.*A;
    return
elseif ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'intop' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end
if (~B.chkval && ~isvalid(B))
    error("The functional is not properly declared.")
end

% Allow for multiplication of scalar with matrix
Adim = size(A);     Bdim = size(B);
if all(Adim==1)
    if isa(A,'double')
        C = B;
        C.params = cellfun(@(a) a*A,B.params,'UniformOutput',false);
    else
        C = intop(cellfun(@(a) A*a,kernels(B),'UniformOutput',false),B.intmat,B.vars,B.dom);
    end
    return
end

% Make sure the row dimensions of the objects match
if Adim(1)~=Bdim(1) && Adim(1)~=1 && Bdim(1)~=1
    error("Row dimensions of objects to multiply must match.")
end
% Make sure the column dimensions of the objects match
if Adim(2)~=Bdim(2) && Adim(2)~=1 && Bdim(2)~=1
    error("Column dimensions of objects to multiply must match.")
end

% Finally, multiplication with the functional amounts to multiplication
% with the parameters
if ~isa(A,'double')
    % Polynomial or dpvar multiplier: multiply the kernels, and express them in
    % terms of the same monomials and decision variables again
    C = intop(cellfun(@(a) A.*a,kernels(B),'UniformOutput',false),B.intmat,B.vars,B.dom);
    return
end
% Constant multiplier: first let the functional take the dimensions of the
% product, as a singleton dimension is repeated in the elementwise product
C = B;
nd = numel(B.Zd);
Cdim = max(Adim,Bdim);
if Bdim(1)==1 && Cdim(1)>1
    C.params = cellfun(@(a) lmult_params(ones(Cdim(1),1),a,nd),C.params,'UniformOutput',false);
end
if Bdim(2)==1 && Cdim(2)>1
    C.params = cellfun(@(a) repmat(a,1,Cdim(2)),C.params,'UniformOutput',false);
end
% Multiply the nonzeros by the entries of A
A = A.*ones(Cdim);
[~,~,nZ] = Zdims(B.Z.in);
C.params = cellfun(@(a) times_param(a,A,nd,nZ,Cdim(1)),C.params,'UniformOutput',false);

end


function P = times_param(P,A,nd,nZ,m)
% Multiply the coefficient matrix P by the elements of the matrix A
[i,j,v] = find(P);
i = i(:);   j = j(:);   v = v(:);       % row vectors if P has one row
row = ceil(i/(nd+1));
col = ceil(j/nZ);
Av = A(:);                              % column, also if A is a row vector
v = v.*Av(row + (col-1)*m);
P = sparse(i,j,v,size(P,1),size(P,2));
end
