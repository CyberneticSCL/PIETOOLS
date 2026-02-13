function C = mtimes(A,B)
% C = MTIMES(A,B) returns the 'intop' object C representing the product of
% the multiplier operator defined by A with the functional defined by B.
%
% INPUTS
% - A:      m x p 'double', 'polynomial', or 'dpvar' object representing 
%           some multiplier operator to multiply with B;
% - B:      p x n 'intop' object acting on some degree-d distributed 
%           monomial;
%
% OUTPUTS
% - C:      m x n 'intop' object representing the product of the
%           multiplier A with the functional defined by B;
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
% DJ, 02/03/2026: Initial coding

% Check that multiplication is supported
if isa(A,'intop')
    error("Left-multiplication with 'intop' objects is not supported.")
elseif ~isa(A,'dobule') && ~isa(A,'polynomial') && ~isa(A,'dpvar')
    error("'intop' objects can only be multiplied with objects of type 'double', 'polynomial', or 'dpvar'.")
end

% Check that the inner dimensions of the objects match
Adim = size(A);
Bdim = size(B);
if Adim(2)~=Bdim(1) && ~all(Adim==1)
    error("Inner dimensions of the objects to multiply must match.")
end

% Multiplication with the functional amounts to multiplication with the
% parameters
C = B;
C.params = A*B.params;

end