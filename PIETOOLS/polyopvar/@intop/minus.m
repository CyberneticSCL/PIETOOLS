function C = minus(A,B)
% C = MINUS(A,B) returns the 'intop' object C representing the
% difference betweeen the functionals defined by the 'intop' objects A and 
% B
%
% INPUTS
% - A:      m x n 'intop' object acting on a degree-d distributed
%           monomial;
% - B:      m x n 'intop' object acting on the same degree-d distributed 
%           monomial as A;
%
% OUTPUTS
% - C:      m x n 'intop' object representing the difference between the 
%           functionals defined by A and B;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - minus
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

C = plus(A,uminus(B));

end