function PDE = subsasgn(PDE,prop,val)
% SUBSASGN allows properties and elements of the "pde_struct" object PDE to
% be adjusted or set. The function is automatically called when calling
% e.g. PDE(i,j) = val or PDE.obj = val for a "pde_struct" PDE.
%
% INPUTS:
% - PDE:    A pde_struct object of which to adjust a property/element.
% - prop:   A 1xq struct specifying which property/element of the input PDE
%           structure to adjust.
% - val:    Value to assign to the specified property/element of the PDE.
%
% OUTPUTS:
% - PDE:    A pde_struct object almost identical to the input system, but
%           with the element/prop specified by "prop" now set equal to
%           "val".
%
% EXAMPLE: calling PDE.x{1}.dom(2) = b, we have
%   PDE_in = PDE;
%   prop(1).type = '.',     prop(1).subs = 'x';
%   prop(2).type = '{}',    prop(2).subs = {[1]};
%   prop(3).type = '.',     prop(3).subs = 'dom';
%   prop(4).type = '()',    prop(4).subs = {[2]};
%   val = b;
% The returned object PDE will be the same as the input object, only with 
% PDE.x{1}.dom(2) = b.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2022 PIETOOLS Team
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
% Initial coding DJ - 10/17/2022
% DJ, 01/03/2025: Account for fact that PDE variables are now represented
%                   by a free term (see also the update to pde_var);
% DJ, 01/04/2025: Allow order of differentiability of state to be set;
%

% At this point, we just use the built-in subsasgn function. Additional
% features such as '()' type subsasgn may be introduced in later updates.
if strcmp(prop(1).type,'()') || strcmp(prop(1).type,'{}')
    error(['''()'' and ''{}'' type subsasgn are currently not supported for "pde_struct" class objects.'])
end

% If the PDE structure corresponds to a single variable (state, input, or
% output), allow the variables, domain, and size to be set.
[is_var,obj] = is_pde_var(PDE);
if is_var && strcmp(prop(1).type,'.') && ismember(prop(1).subs,{'size';'vars';'var';'dom';'diff'})
    % Allow the size of the object to be set.
    if strcmp(prop(1).subs,'size')
        if ~isnumeric(val) || numel(val)~=1 || val<=0 || round(val)~=val
            error("The desired size of the PDE variable should be specified as a scalar positive integer.")
        else
            PDE.(obj){1}.size = val;
            PDE.([obj,'_tab'])(1,2) = val;
            PDE.free{1}.size = val;                                         % DJ, 01/03/2025
            PDE.free{1}.term{1}.C = eye(val);
        end
    elseif strcmp(prop(1).subs,'vars') || strcmp(prop(1).subs,'var')
    % Allow the spatial variables on which the object depends to be set.
        if ischar(val) && size(val,1)==1
            % Convert variable name to polynomial.
            vars = polynomial({val});
        elseif iscellstr(val)
            % Convert variable names to polynomials.
            vars = polynomial(val);
        elseif isempty(val)
            % Set empty list of variables.
            vars = polynomial(zeros(0,1));
        elseif ~isempty(val) && ~ispvar(val)
            error('Spatial variables for the desired object should be specified as an object of type ''polynomial'' or as a cell of variable names.' )
        else
            vars = val;
        end
        if ~isempty(val) && ~any(size(vars)==1)
            error('Spatial variables for the desired object should be specified as array of size nx1.')
        else
            vars = vars(:);
        end
        % Remove any temporal variables, if specified.
        if ismember('t',vars.varname)
            vars = vars(~isequal(vars,polynomial({'t'})));
        end
        PDE.(obj){1}.vars = vars;
        PDE.free{1}.vars = vars;                                            % DJ, 01/03/2025
        PDE.free{1}.term{1}.loc = vars';
        PDE.free{1}.term{1}.D = zeros(1,size(vars,1));
        PDE.free{1}.term{1}.I = cell(size(vars,1),1);
        % Also set a domain for the variables.
        PDE.(obj){1}.dom = [zeros(size(vars,1),1),ones(size(vars,1),1)];
    elseif strcmp(prop(1).subs,'dom')
    % Allow the domain of the variables to be set.
        % If variables are specified, the number of domains should match
        % the number of variables.
        if isfield(PDE.(obj){1},'vars')
            nvars = size(PDE.(obj){1}.vars,1);
        else
            nvars = nan;
        end
        dom = val;
        if ~isa(dom,'double')
            error('The domain of the spatial variables should be specified as nx2 array for n variables.')
        elseif size(dom,2)~=2
            error('The domain of the spatial variables should be specified as nx2 array for n variables.')
        elseif any(dom(:,1)>=dom(:,2))
            error('The specified domain makes no sense: lower boundary exceeds upper boundary!')
        elseif ~isnan(nvars)
            if size(dom,1)==1
                % Assume same domain for all vars
                dom = repmat(dom,nvars(1),1);
            elseif size(dom,1)~=nvars
                error('The number of specified domains should match the number of specified variables.')
            end
        end
        PDE.(obj){1}.dom = dom;
    elseif strcmp(prop(1).subs,'diff')                                      % DJ, 01/04/2025
    % Allow the order of differentiability of a state wrt its spatial
    % variables to be set.
        if ~isnumeric(val) || any(val<0) || any(round(val)~=val)
            error("Order of differentiability should be specified as nx1 array of nonnegative integers.")
        elseif ~strcmp(obj,'x')
            error("Setting order of differentiability of input or output signals is not supported.")
        end
        % Make sure an order of differentiability is specified for each
        % variable.
        if isfield(PDE.(obj){1},'vars') && numel(val)~=size(PDE.(obj){1}.vars,1)
            if isscalar(val)
                % Assume same order of differentiability wrt all vars
                val = val*ones(1,size(PDE.(obj){1}.vars,1));
            else
                error("For a state in n spatial variables, the order of differentiability should be specified as nx1 array.")
            end
        end
        val = val(:)';
        PDE.(obj){1}.diff = val;
    end
    return
else
    PDE = builtin('subsasgn',PDE,prop,val);
end


% If any property of the PDE has been adjusted, the PDE is no longer
% initialized (in case it was before).
if PDE.is_initialized && ~strcmp(prop(1).subs,'is_initialized')
   PDE.is_initialized = false;
end

end