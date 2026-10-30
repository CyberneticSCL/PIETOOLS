function [vars,lims,is_var] = parse_int_inputs(vars,dom,vars_reg,dom_reg)
% [VARS,LIMS,IS_VAR] = PARSE_INT_INPUTS(VARS,DOM,VARS_REG,DOM_REG) checks
% and normalizes the variables and limits of integration specified for the
% 'int' functions of the 'sopvar' and 'copvar' classes, and 'int2sopvar'.
%
% INPUTS
% - vars:       N x 1 or 1 x N 'cellstr', char, string or 'polynomial'
%               (pvar) array specifying the unique names of the variables
%               to integrate over;
% - dom:        (optional) N x 2 'double', 'polynomial' or 'cell' array,
%               dom(i,:) = [L,U] the lower and upper limit of integration
%               in vars(i). A limit is a number, or the variable vars(i)
%               itself (a pvar, or a char/string in a cell). Defaults to
%               the full domain of each variable, as given by vars_reg and
%               dom_reg;
% - vars_reg:   1 x nv 'cellstr' of the variables with a known domain;
% - dom_reg:    nv x 2 'double', dom_reg(k,:) = [a,b] the domain of
%               vars_reg{k};
%
% OUTPUTS
% - vars:       1 x N 'cellstr' of the names of the variables;
% - lims:       N x 2 'double' of the numeric limits, NaN where the limit
%               is the variable itself;
% - is_var:     N x 2 'logical', true where the limit is the variable
%               itself;
%
% NOTES
% - Whether a limit is admissible depends on the operator. It is checked in
%   'int2sopvar', and in the 'int' functions for the variables that are not
%   variables of the operator;
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PIETOOLS - parse_int_inputs
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
% DJ (CC), 10/29/2026: Initial coding. Shared by @sopvar/int, @copvar/int
%                      and int2sopvar, so that they accept the same inputs.
%                      Unlike parse_diff_inputs, the names are not sorted:
%                      row i of dom belongs to vars{i}.

% Variables: normalize to a 1 x N cellstr, so one list serves all types
if isa(vars,'polynomial')
    if ~ispvar(vars)
        error("Integration with respect to polynomial combinations of variables is not supported.")
    end
    vars = pvar2varname(vars);
elseif ischar(vars)
    vars = {vars};
elseif isstring(vars)
    vars = cellstr(vars);
elseif ~iscellstr(vars)
    error("Variable names should be specified as 'cellstr' or 'polynomial' array.")
end
vars = reshape(vars,1,[]);      N = numel(vars);
% Integrating twice in one variable is not a single operator in that variable
if numel(unique(vars))<N
    error("Variables of integration must be unique.")
end

% Default: the full domain, so every limit is a number from the registry
if nargin<=1 || isempty(dom)
    [is_reg,loc] = ismember(vars,vars_reg);
    if ~all(is_reg)
        error("The domain of variable '%s' is unknown, so it must be specified.",vars{find(~is_reg,1)})
    end
    lims = dom_reg(loc,:);      is_var = false(N,2);
    return
end

% Limits: normalize to an N x 2 cell, so each entry is read the same way
if isa(dom,'polynomial') || isnumeric(dom)
    if ~isequal(size(dom),[N,2])
        error("Domain must be specified as N x 2 array, with a row for each of the N variables.")
    end
    if isnumeric(dom)
        dom = num2cell(dom);
    else
        dom = arrayfun(@(k) dom(k),reshape(1:2*N,N,2),'UniformOutput',false);
    end
elseif ~iscell(dom) || ~isequal(size(dom),[N,2])
    error("Domain must be specified as N x 2 'double', 'polynomial' or 'cell' array, with a row for each of the N variables.")
end
% Read each limit: a number, or the variable itself, as pvar or by name. A
% limit in any other variable would make the output depend on that variable
% through the limit, which is not an operator on the same spaces
lims = zeros(N,2);      is_var = false(N,2);
for i=1:N
    for j=1:2
        lim_ij = dom{i,j};
        vname_lim = '';
        if isa(lim_ij,'polynomial')
            if ispvar(lim_ij)
                vname_lim = pvar2varname(lim_ij);       vname_lim = vname_lim{1};
            elseif isempty(lim_ij.varname)
                lim_ij = double(lim_ij);
            else
                error("Limits of integration in '%s' must be numbers or '%s' itself.",vars{i},vars{i})
            end
        elseif ischar(lim_ij) || isstring(lim_ij)
            vname_lim = char(lim_ij);
        end
        if ~isempty(vname_lim)
            if ~strcmp(vname_lim,vars{i})
                error("Limit of integration '%s' not allowed in '%s': limits must be numbers or '%s' itself.",vname_lim,vars{i},vars{i})
            end
            lims(i,j) = NaN;        is_var(i,j) = true;
            continue
        end
        % An infinite limit would make the constant factor (U-L) infinite
        if ~isnumeric(lim_ij) || ~isscalar(lim_ij) || ~isreal(lim_ij) || ~isfinite(lim_ij)
            error("Limits of integration in '%s' must be finite real numbers or '%s' itself.",vars{i},vars{i})
        end
        lims(i,j) = lim_ij;
    end
end

end
