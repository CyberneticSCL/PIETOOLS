function C = kron(A,B)
% C = KRON(A,B) takes the Kronecker product of two matrix-valued tensor-PI
% operators A and B.
%
% INPUTS
% - A:  m x n 'tensopvar' object (or 'sopvar', converted to 'tensopvar');
% - B:  p x q 'tensopvar' object (or 'sopvar', converted to 'tensopvar');
%
% OUTPUTS
% - C:  m*p x n*q 'tensopvar' object, representing the Kronecker product of
%       A and B;
%
% NOTES
% - This function does not compute the tensor product of the operators,
%   use `otimes` for that instead.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kron
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

% Check the inputs
if isa(A,'sopvar')
    A = sopvar2tensopvar(A);
elseif ~isa(A,'tensopvar') && ~isnumeric(A)
    error("Inputs must be numeric or of type 'tensopvar'.")
end
if isa(A,'tensopvar') && ~A.chkval
    A = initialize(A);
end
if isa(B,'sopvar')
    B = sopvar2tensopvar(B);
elseif ~isa(B,'tensopvar') && ~isnumeric(B)
    error("Inputs must be numeric or of type 'tensopvar'.")
end
if isa(B,'tensopvar') && ~B.chkval
    B = initialize(B);
end

% Compute the kronecker product using a rudimentary approach
C = otimes(repelem(A,size(B)),repmat(B,size(A)),[1,1]);

end