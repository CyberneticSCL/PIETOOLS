function [logval,varname,dims,vars,dom] = ispvar(p)
% [LOGVAL,VARNAME,DIMS,VARS,DOM] = ISPVAR(P) checks whether the
% 'polyopvar' object P is an array of individual (state) variables,
% returning an array of the state variables and the spatial variables that
% they depend on.
%
% INPUTS
% - p:          mtot x 1 'polyopvar' object;
%
% OUTPUTS
% - logval:     boolean indicating whether (true) or not (false) each of
%               the elements of p corresponds to an individual state
%               variable (acted upon by an identity operator);
% - varname:    m x 1 cellstr, where m = size(p.C.ops,1) is the number of
%               blocks of rows of p.C, specifying the names of the state
%               variables in each row of p, if logval==true. Otherwise,
%               varname = cell(0,1). Elements may not be unique: the same
%               state variable may appear in different rows in p;
% - dims:       m x 1 array of integers, specifying the vector size of each
%               of the state variables listed in "varname". If
%               logval==false, dims = zeros(0,1);
% - vars:       m x 1 cell with each element vars{i} a 1 x ni cellstr
%               specifying the names of the spatial variables on which
%               state varname{i} depends. If logval==false, vars =
%               cell(0,1);
% - dom:        m x 1 cell with each element a ni x 2 array specifying for
%               each of the variables in vars{i} the interval on which
%               these variables are defined. If logval==false, dom =
%               cell(0,1);
%
% NOTES
% - Polynomial arrays with multiple columns are not supported.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - ispvar
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

% Make sure the polynomial is properly declared
if ~p.chkval
    p = initialize(p);
end

% Assume p is not pvar until we have checked
logval = false;
varname = cell(0,1);
dims = zeros(0,1);
vars = cell(0,1);
dom = cell(0,1);

% Check if the degrees of the monomials are all equal to 1
is_nz = p.Z>0;
if any(sum(is_nz,2)~=1)
    return
end
Z = p.Z';           % Z = p.Z(is_nz) may mess up the order
Z = Z(is_nz');
Z = Z';

% Check that each of the rows of p is defined by an identity operator
% acting on some state
% Block dimensions of C: size(p.C) is the matrix dimension, the sum of the sizes of the blocks
[nr,nc] = size(p.C.ops);
Zidcs = zeros(nr,1);
for i=1:nr  
    % Find nonzero column of C
    for j=1:nc
        if isempty(p.C{i,j})
            % C{i,j} is zero
            continue
        end
        % Check if C{i,j} is a zero operator
        [Cij,is_zero_ij] = clean_op(p.C{i,j});
        if is_zero_ij
            continue
        end
        % Check if Cij is an identity operator
        is_id_ij = is_id_op(Cij);
        if ~is_id_ij || Zidcs(i)~=0
            % Row i involves a non-identity operator, or multiple state
            % variables
            return
        else
            Zidcs(i) = j;
        end
    end
    if Zidcs(i)==0
        return
    end
end

% At this point, we can conlude that all rows of p involve a single state
% --> identify this state
logval = true;
var_idcs = Z(Zidcs);
varname = reshape(p.varname(var_idcs(:)),[],1);
dims = p.dims.in(var_idcs);
[var_idcs_unique,~,new2old_idcs] = unique(var_idcs(:));
vars = cell(numel(var_idcs_unique),1);
dom = cell(numel(var_idcs_unique),1);
% Also determine on what variables each state depends
for i=1:numel(var_idcs_unique)
    deps_i = p.depmat.in(var_idcs_unique(i),:);
    vars{i} = p.vars.in(deps_i);
    dom{i} = p.dom.in(deps_i,:);
end
vars = vars(new2old_idcs);
dom = dom(new2old_idcs);

end