function B = ctranspose(A)
% B = CTRANSPOSE(A) would return the adjoint A' of the 'intop' object
% A. The adjoint of these operators is currently not supported, so an
% error is thrown.
%
% INPUTS
% - A:      'intop' object;
%
% OUTPUTS
% - B:      none, an error is thrown;
%
% NOTES
% - The adjoint of a distributed polynomial U, with U'*V = innerprod(U,V),
%   is supported for 'polyopvar' objects, see 'polyopvar/ctranspose'.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ctranspose
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
% DJ, 11/07/2026: Initial coding;

error("intop:adjoint_not_supported", ...
    "The adjoint of 'intop' objects is currently not supported.")

end
