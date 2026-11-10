function f = remove_unused_states(f)
% F = REMOVE_UNUSED_STATES(F) takes a 'polyopvar' object F and removes the
% state variables that do not appear in any of its monomials, remapping
% F.Z to the remaining variables.
%
% INPUTS
% - f:      'polyopvar' object representing a distributed polynomial
%           P(x) = C*Z(x);
%
% OUTPUTS
% - f:      'polyopvar' object representing the same polynomial, in
%           terms of only those state variables (entries of f.varname)
%           that appear in some monomial (row of f.Z);
%
% NOTES
% - Only f.Z, f.varname, f.dims.in and f.depmat.in are modified: the
%   coefficients are left as is, in particular no tolerance is applied
%   (see clean_poly for that).
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - remove_unused_states
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

% Determine which state variables appear in a monomial
p = numel(f.varname);
[ridcs,cidcs,vals] = find(f.Z);
state_used = false(p,1);
state_used(unique(vals)) = true;
if all(state_used)
    return
end

% Remove the unused variables, and remap the indices in Z
remap = zeros(p,1);
remap(state_used) = 1:nnz(state_used);
f.Z = sparse(ridcs,cidcs,remap(vals),size(f.Z,1),size(f.Z,2));
f.varname = reshape(f.varname(state_used),[],1);
f.dims.in = f.dims.in(state_used);
f.depmat.in = f.depmat.in(state_used,:);
f.chkval = true;

end
