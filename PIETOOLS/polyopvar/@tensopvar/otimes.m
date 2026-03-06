function C = otimes(A,B)
% C = OTIMES(A,B) returns the tensor product of coefficient operators A and
% B, so that (A*Zd(x))(s).*(B*Zd(y))(s) = (C*(Zd(x) o Zd(y)))(s)
%
% INPUTS
% - A:  m x n1 'tensopvar' object
% - B:  m x n2 'tensopvar' object
%
% OUTPUTS
% - C:  m x n1*n2 'tensopvar' object representing the tensor product of A
%       and B
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
% DJ, 03/01/2026: Initial coding


if size(A,1)~=size(B,1)
    error("Row-dimension of operators to multiply must match.")
end

% Allow for multiplication with empty operator
d1 = numel(A.ops);
d2 = numel(B.ops);
if d1==0
    C = B;
    return
elseif d2==0
    C = A;
    return
end

% Construct the operator acting on each product of distributed monomials
Ccell = cell(1,d1*d2);
for k=1:numel(Ccell)
    % Determine which term in A and B corresponds to term k in C
    [k1,k2] = ind2sub([d1,d2],k);
    Ak = A.ops{k1};
    Bk = B.ops{k2};

    % Take the tensor product of the operators Ak and Bk
    if isa(Ak,'cell') && isa(Bk,'cell')
        % Multiplication of multiple factors in polynomial vector field
        % --> just concatenate the operators
        Ck = [repmat(Ak,size(Bk,1),1),repelem(Bk,size(Ak,1),1)];
    elseif isa(Ak,'cell')
        % Multiplication of multiple factors with single factor
        if isa(Bk,'intop')
            error("Multiplication of polynomial functional with polynomial vector field is not supported.")
        elseif ~isa(Bk,'nopvar')
            Bk = nopvar(Bk);
        end
        Ck = [Ak,repmat({Bk},size(Ak,1),1)];
    elseif isa(Bk,'cell')
        % Multiplication of single factor with multiple factors
        if isa(Ak,'intop')
            error("Multiplication of polynomial functional with polynomial vector field is not supported.")
        elseif ~isa(Ak,'nopvar')
            Ak = nopvar(Ak);
        end
        Ck = [repmat({Ak},size(Bk,1),1),Bk];
    else
        if isa(Ak,'intop') || isa(Bk,'intop')
            Ck = otimes(Ak,Bk);
        elseif isa(Ak,'nopvar') && isa(Bk,'nopvar')
            % Tensor product of two PI operators
            Ck = {Ak,Bk};
        else
            % Multiplication with degree-0 term
            Ck = Ak.*Bk;            
        end
    end
    Ccell{k} = Ck;
end

% Declare the full coefficient operator
C = tensopvar();
C.ops = Ccell;

end