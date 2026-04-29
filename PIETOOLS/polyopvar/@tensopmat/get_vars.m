function [vars,dom] = get_vars(Cop)
% [VARS,DOM] = GET_VARS(COP) returns the spatial and dummy variables in
% terms of which the operator COP is defined
%
% INPUTS
% - Cop:    'tensopmat' object representing a matrix of tensor-PI operators
%           mapping between different function spaces;
%
% OUTPUTS
% - vars:   N x 2 'pvar' array specifying the spatial variables (first
%           column) and dummy variables (second column) in terms of which
%           the operator is defined;
% - dom:    N x 2 array specifying for each of the spatial variables the
%           interval [a,b] on which this variable is defined;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_vars
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
% DJ, 04/22/2026: Initial coding

% Initialize empty arrays
varname1 = {};
varname2 = {};
dom = zeros(0,2);
for j=1:numel(Cop.ops)
    % Loop over all tensor-PI operators in the structure
    Cop_j = Cop.ops{j};
    if isempty(Cop_j)
        continue
    end
    % Extract the properties of the tensor-PI operator in block j
    vars_j = Cop_j.vars;
    var1_j = pvar2varname(vars_j(:,1));
    var2_j = pvar2varname(vars_j(:,2));
    dom_j = Cop_j.dom;
    % Add to our full lists of properties
    varname1 = [varname1; var1_j];
    varname2 = [varname2; var2_j];
    dom = [dom;  dom_j];
    % Get rid of duplicates
    [varname1,idcs] = unique(varname1,'stable');
    varname2 = varname2(idcs);
    dom = dom(idcs,:);
end
vars = polynomial([varname1,varname2]);

end