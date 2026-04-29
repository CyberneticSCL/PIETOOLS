function [depmat1,depmat2] = get_deps(Cop)
% [DEPMAT1,DEPMAT2] = GET_DEPS(COP) returns the input and output dependency
% arrays of COP in terms of the variables in which this COP is defined.
%
% INPUTS
% - Cop:    'tensopmat' object representing a matrix of tensor-PI operators
%           mapping between different function spaces;
%
% OUTPUTS
% - depmat1:    nr x N array where nr is the number of rows of operators in
%               Cop, and N the number of spatial variables. Then
%               depmat1(i,j) indicates for row i of the output of Cop to
%               what degree it depends on spatial variable j (where degree
%               refers to the order of tensor product along this spatial
%               variable);
% - depmat1:    nc x N array where nc is the number of columns of operators
%               in Cop, and N the number of spatial variables. Then
%               depmat2(i,j) indicates for column i of the input of Cop to
%               what degree it depends on spatial variable j (where degree
%               refers to the order of tensor product along this spatial
%               variable);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - get_deps
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

% Extract spatial variables used in definition of the operator
varname = pvar2varname(Cop.vars(:,1));
nvars = numel(varname);
% Initialize empty dependency arrays
[nr,nc] = size(Cop.ops);
depmat1 = -1*ones(nr,nvars);
depmat2 = -1*ones(nc,nvars);
for lidx=1:numel(Cop.ops)
    % Loop over all tensor-PI operators in the structure
    if isempty(Cop.ops{lidx})
        continue
    end
    [ridx,cidx] = ind2sub(size(Cop.ops),lidx);
    var1_j = pvar2varname(Cop.ops{lidx}.vars(:,1));
    deps = Cop.ops{lidx}.dep;
    % For each of the variables appearing in the tensor-PI operator,
    % set the dependence on this variable in the tensopmat structure
    for k=1:numel(var1_j)
        var_idx = find(ismember(varname,var1_j{k}),1);
        if depmat1(ridx,var_idx)<0
            depmat1(ridx,var_idx) = deps(1,k);
        elseif depmat1(ridx,var_idx)~=deps(1,k)
            depmat1(ridx,var_idx) = nan;
        end
        if depmat2(cidx,var_idx)<0
            depmat2(cidx,var_idx) = deps(2,k);
        elseif depmat2(cidx,var_idx)~=deps(2,k)
            depmat2(cidx,var_idx) = nan;
        end
    end    
end
% Assume empty operators map from R to R (no spatial variables)
depmat1(depmat1<0) = 0;
depmat2(depmat2<0) = 0;

end