function At = transpose(A)
% At = TRANSPOSE(A) returns the adjoint of the 'sdopvar' operator A,
% A.' = A'; it simply calls 'ctranspose', which swaps the output and input
% variables and the integral kernels as described there.
%
% INPUTS
% - A:  'sdopvar' object representing a PI operator L_2^p[S1,S3] to
%       L_2^q[S2,S3];
%
% OUTPUTS
% - At: 'sdopvar' object representing the adjoint operator, mapping
%       L_2^q[S2,S3] to L_2^p[S1,S3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - transpose
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
% Initial coding MP, SS, DJ, AT, DB - 10/24/2026
%

At = ctranspose(A);
end
