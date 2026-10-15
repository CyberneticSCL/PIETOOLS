function nkey = factor_key(Z,nops,is_empty)
% NKEY = FACTOR_KEY(Z,NOPS,IS_EMPTY) returns the numbers of factors to use,
% together with the rows of Z, to identify equal monomials.
%
% INPUTS
% - Z:          nZ x d array, the monomials of a 'polyopvar', or of
%               several, stacked;
% - nops:       nZ x 1 array, the number of factors in the operators
%               acting on each monomial, see num_factors;
% - is_empty:   nZ x 1 logical array, true if all operators acting on the
%               monomial are empty (zero), see num_factors;
%
% OUTPUTS
% - nkey:       nZ x 1 array of integers. Monomials with the same Z but
%               different numbers of factors are different monomials,
%               e.g. x o c (c a multiplier from R, Z=[1,0], two factors)
%               and x (Z=[1], one factor), which are the same row of Z
%               once padded to the same width. A monomial that is not
%               acted on by any operator has no factors to distinguish,
%               so it takes the number of factors of a monomial with the
%               same Z that is acted on by operators (the first), and so
%               can be merged with it;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - factor_key
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

nkey = nops;
if any(is_empty) && any(~is_empty)
    idcs_e = find(is_empty);
    idcs_n = find(~is_empty);
    [is_found,loc] = ismember(full(Z(idcs_e,:)),full(Z(idcs_n,:)),'rows');
    nkey(idcs_e(is_found)) = nops(idcs_n(loc(is_found)));
end

end
