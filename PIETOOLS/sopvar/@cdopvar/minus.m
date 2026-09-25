function P = minus(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% P = MINUS(A,B) or A-B subtracts containers blockwise.
%
% INPUTS
% - A, B:   'cdopvar' or 'copvar' objects on the same grid, spaces and
%           dimensions, as 'plus' requires; at least one is a 'cdopvar'
%           for this method to be reached;
%
% OUTPUTS
% - P:      'cdopvar' object equal to A - B.
%
% NOTES
% Written as A + (-B), as @sdopvar/minus is, so that the promotion of a
% 'copvar' operand and the reconciliation of decision variable lists live
% only in @cdopvar/plus. Negating B leaves its list unchanged.
%
% Cost: that of 'uminus' on B plus that of 'plus'.
%
% See also UMINUS, PLUS, CDOPVAR.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - minus(cdopvar)
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

P = plus(A,uminus(B));

end
