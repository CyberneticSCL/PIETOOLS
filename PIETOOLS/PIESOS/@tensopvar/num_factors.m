function n = num_factors(Cop)
% N = NUM_FACTORS(COP) returns the number of factors in the 'tensopvar'
% object COP, counting also the factors of nested 'tensopvar' objects.
%
% INPUTS
% - Cop:    'tensopvar' object representing a tensor-PI operator
%               C = sum_{i=1}^{K} Ri1 \otimes ... \otimes Rid;
%
% OUTPUTS
% - n:      number of factors in the first term, d if all factors are
%           'sopvar' or 'intop' objects;
%
% NOTES
% - The number of factors is that of the first term; 'initialize' checks
%   that all terms have the same structure.
% - This is a method of the class, so that Cop.ops is read directly as a
%   property, not through the 'subsref' of the class (see
%   @polyopvar/num_factors.m).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - num_factors (tensopvar)
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

fctrs = Cop.ops(1,:);
n = numel(fctrs);
for j=find(cellfun('isclass',fctrs,'tensopvar'))
    n = n - 1 + num_factors(fctrs{j});
end

end
