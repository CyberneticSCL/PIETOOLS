function Cop = minus(Aop,Bop)
% COP = MINUS(AOP,BOP) returns the nopvar object representing the
% difference of the PI operators defined by AOP and BOP.
%
% INPUTS
% - Aop:    m x n 'nopvar' object
% - Bop:    m x n 'nopvar' object
%
% OUTPUS
% - Cop:    m x n 'nopvar' object representing the sum of the operators
%           defined by Aop and Bop
%
% NOTES
% The operators defined by Aop and Bop must act on functions on the same
% spatial domain, i.e. Aop.dom = Bop.dom.

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
% DJ, 01/15/2026: Initial coding

Cop = plus(Aop,uminus(Bop));

end