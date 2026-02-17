function PDE_out = times(fctr1,fctr2)
% PDE_OUT = TIMES(FCTR1,FCTR2) declares a new PDE term corresponding to
% the elementwise product of the factors FCTR1 and FCTR2.
%
% INPUT
% - fctr1:      'pde_struct' object, 'polynomial', or 'double',
%               representing a vector of coefficients or a set of terms of 
%               a PDE;
% - fctr2:      'pde_struct' object representing a set of PDE terms to be
%               multiplied by fctr1. The number of rows of terms, as well
%               as the vector dimensions of the terms in each row, must
%               match those of fctr1 if fctr1 is a 'pde_struct' object as
%               well.
%
% OUTPUT
% - PDE_out:    'pde_struct' object representing one or several terms in a
%               PDE, corresponding to the elementwise product fctr1.*fctr2;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026  PIETOOLS Team
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% DJ, 02/17/2026: Initial coding;

% % % Process the inputs

% % Check that the first input makes sense.
if ~isa(fctr2,'pde_struct')
    error("Right-multiplication of 'pde_struct' objects with non-'pde_struct' objects is not supported.")
end
if isa(fctr1,'pde_struct')
    % Use separate function for multiplication of PDE terms
    PDE_out = times_pdes(fctr1,fctr2);
    return
elseif ~isa(fctr1,'double') && ~isa(fctr1,'polynomial') && ~isa(fctr1,'opvar')
    error("The first factor in the product should be of type 'pde_struct', 'double' or 'polynomial'.")
elseif ~isa(fctr1,'opvar')
    fctr1 = polynomial(fctr1);
end

% For now, only allow multiplication of scalar coefficient with PDE
if all(size(fctr1)==1)
    PDE_out = mtimes(fctr1,fctr2);
else
    error("Elementwise multiplication of coefficients with PDE terms is currently not supported.")
end

end