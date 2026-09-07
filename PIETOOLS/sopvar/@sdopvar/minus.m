function Cop = minus(Aop,Bop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cop = minus(Aop,Bop) or Aop-Bop subtracts two sdopvar operators.
%
% INPUT
% Aop,Bop:  'sdopvar' objects of equal dimension, sharing vars and dom. A
%           fixed 'sopvar' operand is accepted, on either side;
%
% OUTPUT
% Cop:      'sdopvar' object representing Aop-Bop;
%
% NOTES:
% Written as Aop + (-Bop), as in @sopvar/minus, so that the reconciliation of
% decision variable lists and monomial bases lives only in @sdopvar/plus.
% The negation is one pass over the nonzeros of the subtrahend and is
% negligible next to that reconciliation.
%
% For support, contact M. Peet, Arizona State University at mpeet@asu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2026 PIETOOLS Team
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
% Initial coding MMP, 09/07/2026

Cop = plus(Aop,uminus(Bop));

end
