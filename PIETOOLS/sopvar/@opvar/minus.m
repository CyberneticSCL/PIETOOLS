function P = minus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = MINUS(A,B) or A-B subtracts two 'opvar' objects blockwise.
%
% INPUTS
% - A, B:   'opvar' objects on the same grid, spaces and dimensions; with a
%           'dopvar' operand the difference is a 'dopvar'. A numeric or
%           'dpvar' scalar c is c*I on a square operator; a 'dpvar' matrix
%           of the total size is a multiplier split by the dimensions;
%
% OUTPUTS
% - P:      'opvar' object equal to A - B.
%
% NOTES
% Computed as A + (-B), so the grids, registries and spaces are checked as
% for a sum.
%
% Cost: one negation of B and one sum.
%
% See also UMINUS, PLUS, OPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - minus(opvar)
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
% Initial coding MMP, 09/25/2026

% A - B as A + (-B)
P = plus(A,uminus(B));

end
