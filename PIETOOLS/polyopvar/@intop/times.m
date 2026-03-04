function C = times(A,B)
% C = TIMES(A,B) returns the 'intop' object C representing the elementwise
% product of the multiplier operator defined by A with the functional 
% defined by B.
%
% INPUTS
% - A:      m x n 'double', 'polynomial', or 'dpvar' object representing 
%           some multiplier operator to multiply with B;
% - B:      m x n 'intop' object acting on some degree-d distributed 
%           monomial;
%
% OUTPUTS
% - C:      m x n 'intop' object representing the elementwise product of
%           the multiplier A with the functional defined by B;
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
if isa(A,'intop') && isa(B,'intop')
    error("Multiplication of functionals is not supported.")
elseif isa(A,'intop')
    % Make sure B is the functional
    C = B.*A;
    return
elseif ~isa(A,'double') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'intop' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end

% Allow for multiplication of scalar with matrix
Adim = size(A);     Bdim = size(B);
if all(Adim==1)    
    C = B;
    C.params = A*B.params;
    return
end

% Make sure the row dimensions of the objects match
if Adim(1)~=Bdim(1)
    if Adim(1)==1
        A = ones(Bdim(1),1).*A;
    elseif Bdim(1)==1
        B.params = ones(Adim(1),1).*B;
    else
        error("Row dimensions of objects to multiply must match.")
    end
end
% Make sure the column dimensions of the objects match
if Adim(2)~=Bdim(2)
    if Adim(2)==1
        A = ones(1,Bdim(2)).*A;
    elseif Bdim(2)==1
        Bparams = zeros(size(B.params,1),0);
        for k=1:size(B.omat,1)
            Bparams = [Bparams,ones(1,Adim(2)).*B.params(:,k)];
        end
        B.params = Bparams;
    else
        error("Column dimensions of objects to multiply must match.")
    end
end

% Finally, multiplication with the functional amounts to multiplication
% with the parameters
C = B;
C.params = A.*B.params;

end