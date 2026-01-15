function Cop = plus(Aop,Bop)
% COP = PLUS(AOP,BOP) returns the 'nopvar' object COP representing the sum
% of the PI operators defined by 'nopvar' objects AOP and BOP
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
% PIETOOLS - plus
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



% Check that the operators can indeed be added
if any(Aop.dim~=Bop.dim)
    error("Dimensions of the operators must match.")
end
if size(Aop.dom,1)~=size(Bop.dom,1) || any(any(Aop.dom~=Bop.dom))
    error("Spatial domains on which operators are defined should match.")
end

% Exclude more complicated addition operations
% --> will need to be included later!
if any(Aop.deg~=Bop.deg)
    error("Addition of operators with different monomial degrees is currently not supported.")
end

% Assuming the same monomial degrees, the sum of the operators is just
% defined by the sum of the coefficients.
Cop = Aop;
for ii=1:numel(Cop.C)
    Cop.C{ii} = Aop.C{ii} + Bop.C{ii};
end

end