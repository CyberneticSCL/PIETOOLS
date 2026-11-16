function varname = pvar2varname(var)
% VARNAME = PVAR2VARNAME(VAR) Takes a pvar object VAR, and returns a cell
% VARNAME specifying the names of the variables appearing in each element
% of VAR
%
% INPUTS
% - var:    m x n 'polynomial' object of which each element is an
%           individual variable;
%
% OUTPUTS
% - varname:    m x n 'cell' with each element specifying the name of the
%               variable in the associated element of "var";
%
% NOTES
% - This function is distinct from using "var.varname", which returns the
%   variable names as a k x 1 cell of unique variable names in alphabetical
%   order;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - nopvar
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
% DJ, 04/14/2026: Initial coding

% Check the input
if isempty(var)
    varname = cell(size(var));
elseif ~ispvar(var)
    error("Input must be a 'polynomial' array of 'pvar' objects.")
end

% Extract components defining the polynomial
Cmat = var.C;
pvarname = var.varname;
degmat = var.degmat;

% Deteremine the variable names in each element
varname = cell(size(var));
for i=1:numel(varname)
    % Determine which monomial appears in this element
    idx = find(Cmat(:,i),1,'first');
    % Determine which variable appears in this monomial
    var_idx = find(degmat(idx,:),1,'first');
    % Set the variable
    varname{i} = pvarname{var_idx};
end

end