function x = pdevar2polyopvar(u)
% X = PDEVAR2POLYOPVAR(U) takes a pde_var object U representing a state 
% variable and returns a polyopvar object X representing the same state
%
% INPUTS
% - u:  1 x 1 'pde_struct' object representing a single state variable;
%
% OUTPUTS
% - x:  1 x 1 'polyopvar' object representing a distributed polynomial
%       variable depending on the same independent (spatial) variables as
%       the input, defined on the same domain.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - pdevar2polyopvar
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
% DJ, 05/08/2026: Initial coding

% Check the input
if ~isa(u,'pde_struct')
    error("Input must be a 'pde_var' object.")
end
[logval,obj] = is_pde_var(u);
if ~logval || ~strcmp(obj,'x')
    error("Input must be a single state variable.")
elseif u.size>1
    error("Conversion of vector-valued states is currently not supported.")
end
% Extract the variables and domain
dom = u.dom;
var1 = u.vars(:,1);

% Declare the polyopvar output
x = polyopvar('x1',var1,dom);

end