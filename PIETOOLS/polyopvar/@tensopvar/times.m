function C = times(A,B)
% C = MTIMES(A,B) returns the 'tensopvar' object C representing the
% elementwise product of the multiplier operator defined by A with the 
% tensor-PI operator defined by B.
%
% INPUTS
% - A:      m x p 'double', 'polynomial', or 'dpvar' object representing 
%           some multiplier operator to multiply with B;
% - B:      p x n 'tensopvar' object acting on some degree-d distributed 
%           monomial;
%
% OUTPUTS
% - C:      m x n 'tensopvar' object representing the product of the
%           multiplier A with the tensor-PI operator defined by B;
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
% DJ, 03/03/2026: Initial coding

% Check that multiplication is supported
if isa(A,'tensopvar') && isa(B,'tensopvar')
    error("Multiplication of tensor-PI operators is not supported.")
elseif isa(A,'tensopvar')
    % Make sure B is the tensor-PI operator
    C = B.*A;
    return
elseif ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'tensopvar' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end


% Allow for multiplication of scalar with matrix
Adim = size(A);     Bdim = size(B);
if all(Adim==1)    
    C = B;
    for j=1:numel(C.ops)
        Cj = C.ops{j};
        if isa(Cj,'cell')
            % Multiply the first factor of each term with A
            for i=1:size(Cj,1)
                Cj{i,1} = A*Cj{i,1};
            end
        else
            % Multiply the operator acting on the jth monomial with A
            Cj = A*Cj;
        end
        C.ops{j} = Cj;
    end
    return
end

% Make sure the row dimensions of the objects match
isdim1 = [false,false];
if Adim(1)~=Bdim(1)
    if Adim(1)==1
        A = ones(Bdim(1),1).*A;
    elseif Bdim(1)==1
        isdim1(1) = true;
    else
        error("Row dimensions of objects to multiply must match.")
    end
end
% Make sure the column dimensions of the objects match
if Adim(2)~=Bdim(2)
    if Adim(2)==1
        A = ones(1,Bdim(2)).*A;
    elseif Bdim(2)==1
        isdim1(2) = true;
    else
        error("Column dimensions of objects to multiply must match.")
    end
end

% Multiply each of the operators stored in B by A
C = B;
for j=1:numel(C.ops)
    Cj = C.ops{j};
    if isa(Cj,'cell')
        % Update dimensions of operator to account for dimensions of B
        if isdim1(1) && isdim1(2)
            for l=1:numel(Cj)
                Cj{l} = ones(Adim(1),Adim(2)).*Cj{l};
            end
        elseif isdim1(1)
            for l=1:numel(Cj)
                Cj{l} = ones(Adim(1),1).*Cj{l};
            end
        elseif isdim1(2)
            for l=1:numel(Cj)
                Cj{l} = ones(1,Adim(2)).*Cj{l};
            end
        end
        % Multiply the first factor of each term with A
        for i=1:size(Cj,1)
            Cj{i,1} = A.*Cj{i,1};
        end
    else
        % Multiply the operator acting on the jth monomial with A
        Cj = A.*Cj;
    end
    C.ops{j} = Cj;
end

end