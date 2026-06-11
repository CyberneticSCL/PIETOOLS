function f_comb = piesos_combine_terms(F)
% f_comb = piesos_combine_terms(F) takes a distributed
% polynomial variable F and canonicalizes repeated distributed monomials 
% using identities like e.g., x(t1)*x(t2) = x(t2)*x(t1).
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial
%           variable;
%
% OUTPUTS
% - prog:   same 'polyopvar' object, but with repeated distributed monomials
%           canonicalized.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - piesos_eq
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
% CRR, 06/09/2026: Initial coding

    f_comb =F;
    for ii = 1:size(F.degmat,1)
        f_ii = F;
        f_ii.degmat = F.degmat(ii,:);
        f_ii.C.ops = F.C.ops(ii);
        f_ii = combine_terms(f_ii); % Make sure to enforce uniqueness of the monomials
        f_comb.C.ops{ii} = f_ii.C.ops{1};
    end
end