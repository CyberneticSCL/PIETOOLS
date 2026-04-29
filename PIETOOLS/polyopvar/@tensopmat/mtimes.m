function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'tensopmat' object C representing the product
% of the multiplier operator defined by A with the tensor-PI operator
% defined by B.
%
% INPUTS
% - A:      m x p 'double', 'polynomial', or 'dpvar' object representing 
%           some multiplier operator to multiply with B;
% - B:      p x n 'tensopmat' object represent a matrix of tensor-PI
%           operators
%
% OUTPUTS
% - C:      m x n 'tensopmat' object representing the product of the
%           multiplier A with the tensor-PI operator defined by B;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - mtimes
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

% Check that multiplication is supported
if isa(A,'nopvar') || isa(A,'tensopvar')
    A = tensopmat(A);
end
if isa(B,'nopvar') || isa(B,'tensopvar')
    B = tensopmat(B);
end
if isa(A,'tensopmat') && isa(B,'tensopmat')
    [A,B] = common_vars(A,B);
    if ~isequal(A.depmat2,B.depmat1)
        error("Input function space of first factor must match output space of second factor.")
    elseif ~isequal(A.dim{2},B.dim{1})
        error("Inner dimensions of operators must match.")
    elseif any(any(A.depmat1>1)) || any(any(A.depmat2>1)) || any(any(B.depmat2>1))
        error("Multiplication of tensor-PI operators is supported only in linear case.")
    end
elseif isa(A,'tensopvar') && all(size(B)==1)
    % Allow right-multiplication with a scalar
    C = B*A;
    return
elseif ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'tensopvar' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end

% Check that the inner dimensions of the objects match
Adim = size(A);
Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Prohibit matrix-multiplication with actual tensor products of PI
% operators for now, since this is not quite trivial
%   (note that  A*(B1*x)(s)*(B2*x)(s) ~= (A*B1*x)(s)*(B2*x)(s))
if ~all(B.dim{1}==1) && any(any(B.depmat1>1))
    error("Matrix-multiplication with tensor-products of PI operators is currently not supported. Use '.*' for elementwise multiplication instead.")
end

if isa(A,'tensopmat')
    % Take the matrix product of the tensor-PI operators
    C = tensopmat();
    C.ops = cell(size(A.ops,1),size(B.ops,2));
    for idx = 1:numel(C.ops)
        [ridx,cidx] = ind2sub(size(C.ops),idx);
        for k=1:size(A.ops,2)
            if ~isempty(A.ops{ridx,k}) && ~isempty(B.ops{k,cidx})
                C.ops{idx} = C.ops{idx} + A.ops{ridx,k}*B.ops{k,cidx};
            end
        end
    end
    % Set the variables and domain of the output
    C.vars = A.vars;
    C.dom = A.dom;
    C.depmat1 = A.depmat1;
    C.depmat2 = B.depmat2;
else
    % Multiply each of the operators stored in B by A
    C = B;
    for i=1:numel(C.ops)
        if ~isempty(B.ops)
            C.ops{i} = A*B.ops{i};
        end
    end
end

end