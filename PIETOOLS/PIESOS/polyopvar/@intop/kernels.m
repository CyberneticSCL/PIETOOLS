function K = kernels(Kop)
% K = KERNELS(KOP) returns the kernels defining the functional KOP, as a
% cell of the same size as KOP.params, with element K{i1,...,iN} the m x n
% kernel K_{i1,...,iN} associated with KOP.params{i1,...,iN}, see the
% documentation of the 'intop' class.
%
% INPUTS
% - Kop:    m x n 'intop' object;
%
% OUTPUTS
% - K:      cell of the same size as Kop.params, with each element an m x n
%           object of type 'dpvar' if the kernel has decision variables, or
%           'polynomial' otherwise, in terms of the dummy variables of the
%           functional, Kop.varnames.in_dum. Each kernel only involves
%           those monomials and those decision variables that appear in it,
%           so that the decision variables of different kernels may differ;
%
% NOTES
% - If Kop has no decision variables, a constant kernel with more than 1e6
%   elements is returned as a sparse array of type 'double' instead.
% - The sparse coefficients Kop.params can be used directly, instead of the
%   'polynomial' and 'dpvar' objects.
% - KERNELS(KOP) is the same as the dependent property KOP.KERNELS, which
%   should be used instead.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - kernels
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
% DJ, 11/04/2026: Return Kop.kernels, see 'get.kernels' in 'intop';

K = Kop.kernels;

end
