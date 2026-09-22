function [deg,report] = degbalance(P,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [DEG,REPORT] = DEGBALANCE(P,OPTS) picks 'possopvar' degrees so that a
% positive operator declared with them carries monomials of similar degree
% to those of P, the n-variate counterpart of '@opvar/degbalance' and
% '@dopvar2d/degbalance'.
%
% The rule, the slot weighting and the limits of what it does are documented
% in 'degbalance_core', which both the 'sopvar' and 'sdopvar' methods share
% because the only class-dependent step is reading a cell's nonzeros.
%
% See also DEGBALANCE_CORE, EQ_OPTS_SOPVAR, POSSOPVAR.
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
% Initial coding MMP, 09/21/2026: thin wrapper so 'degbalance' is reachable
%                  as a method, matching '@opvar', '@dopvar' and
%                  '@dopvar2d'.

if nargin<2
    [deg,report] = degbalance_core(P);
else
    [deg,report] = degbalance_core(P,opts);
end

end
