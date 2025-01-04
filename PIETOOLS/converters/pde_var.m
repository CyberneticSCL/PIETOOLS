function out_var = pde_var(varargin)
% OUT_VAR = PDE_VAR(VARARGIN) generates a "pde_struct" object representing
% a PDE state, input, or output variable. 
% There are two options for using this function:
%
% OPTION 1
% Call e.g.
%   pde_var x1 x2 input w1 w2 output z1 control u1 sense y1 y2
% This will asign to each of the specified variables (x1, x2, w1, etc.) a
% pde_struct object representing a PDE variable. By default, this variable
% will be a state variable, but this can be changed by adding
%   input   - to set subsequent variables to exogenous inputs
%   control - to set subsequent variables to controlled inputs
%   output  - to set subsequent variables to regulated outputs
%   sense   - to set subsequent variables to sensed outputs
%
% OPTION 2
% Use the "standard" input format
%   out_var = pde_var(type,'sz',vars,dom);
%
% INPUT
% - type:   Object of type 'char', specifying what type of PDE variable is
%           desired. Can be one of 'state', 'input', 'control', 'output',
%           or 'sense'. If fewer than 3 input arguments are passed, this
%           defaults to type='state';
% - sz:     Scalar integer specifying the size of the vector-valued
%           variable. Defaults to 1;
% - vars:   nx1 or nx2 object of type 'polynomial' or 'cellstr', specifying
%           the (names of) the n spatial variable on which the declared PDE
%           variable should depend. If specified as nx2 object, the second
%           column should correspond to dummy variables associated to each
%           spatial variable. Defaults to [] (no spatial dependence);
% - dom:    nx2 array of type 'double' specifying for each of the n spatial
%           variables the interval on which it exists. Defaults to
%           repmat([0,1],n,1) for n variables;
%
% OUTPUT
% - out_var:    pde_struct object representing a PDE variable, which can be
%               used to construct a PDE using standard operations of +, *, 
%               subs(), diff(), int(), and ==.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C)2024  M. Peet, S. Shivakumar, D. Jagt
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
% Initial coding DJ - 06/23/2024
% DJ, 12/06/2024: Set default values in case of insufficient arguments;
% DJ, 12/28/2024: Allow length to be omitted as input;
% DJ, 01/03/2025: Declare a "free" term representing the output variable;

% % % Deal with the case the function is used as e.g.
% % %   pde_var x1 x2 x3 input w1 w2
% % % where each input represents a PDE variable to declare
if nargout==0
    % Assume the object must be a state variable, unless specified
    % otherwise
    var_type = 'state';
    for ii=1:nargin
        if ~ischar(varargin{ii})
            error("Each input should be of type 'char', specifying the name of a PDE variables to declare.")
        end
        if strcmpi(varargin{ii},'state') || ...
                strcmpi(varargin{ii},'input') || ...
                strcmpi(varargin{ii},'output') || ...
                strcmpi(varargin{ii},'control') || ...
                strcmpi(varargin{ii},'sense')
            var_type = varargin{ii};
            continue
        end
        % Initialize the object as scalar valued and finitie-dimensional
        obj_ii = pde_var(var_type,1,[],[]);
        assignin('caller',varargin{ii},obj_ii);
    end
    return
elseif nargout>1
    error("At most one output is supported.")
end

% % % Proceed with the case that inputs are specified as
% % %   out_var = pde_var(type,sz,vars,dom)
% Declare default values        % DJ, 12/06/2024
type = 'state';
sz = 1;
vars = [];
dom = [];
if nargin>0 && nargin<=3 && isnumeric(varargin{1})
    % If no type is specified, assume the object corresponds to a 
    % state variable
    if nargin>=1
        sz = varargin{1};
    end
    if nargin>=2
        vars = varargin{2};
    end
    if nargin>=3
        dom = varargin{3};
    end
elseif nargin==2 && isa(varargin{1},'polynomial')                           % DJ, 12/28/2024
    % If no type or length is specified, assume a state of length 1.
    vars = varargin{1};
    dom = varargin{2};
elseif nargin==3 && ischar(varargin{1}) && isa(varargin{2},'polynomial')    % DJ, 12/28/2024
    % If no length is specified, assume a scalar variable.
    type = varargin{1};
    vars = varargin{2};
    dom = varargin{3};
else
    if nargin>=1
        type = varargin{1};
    end
    if nargin>=2
        sz = varargin{2};
    end
    if nargin>=3
        vars = varargin{3};
    end
    if nargin>=4
        dom = varargin{4};
    end
    if nargin>=5
        error('Too many input arguments.')
    end
end

% % % Check that the specified arguments make sense
% % First the desired type of PDE object
if ~ischar(type)
    error("First argument should be of type 'char', corresponding to either 'state', 'input', or 'output'.")
end
% % Then the size of the desired object
if ~isnumeric(sz) || sz<0 || any(sz-round(sz))
    % Declare multiple PDE objects, of lengths sz{k}
    error('The size of the desired variable should be specified as a real, nonnegative, integer value.')
end
% % Then the variables on which the object should depend
if ischar(vars)
    % Convert variable name to polynomial.
    vars = polynomial({vars});
elseif iscellstr(vars)
    % Convert variable names to polynomials.
    vars = polynomial(vars);
elseif isempty(vars)
    % Set empty list of variables.
    vars = polynomial(zeros(0,1));
elseif ~isempty(vars) && ~ispvar(vars)
    error('Spatial variables for the desired object should be specified as an object of type ''polynomial'' or as a cell of variable names.' )
end
if ~isempty(vars) && ~any(size(vars)==1)
    error('Spatial variables for the desired object should be specified as array of size nx1.')
else
    vars = vars(:);
end
nvars = size(vars,1);

% % And finally the domain of each variable.
if isempty(dom)
    % If no domain is specified, assume [0,1]
    dom = [zeros(nvars,1),ones(nvars,1)];
elseif isa(dom,'double')
    if size(dom,2)~=2
        error('The domain of the spatial variables should be specified as nx2 array for n variables.')
    elseif any(dom(:,1)>=dom(:,2))
        error('The specified domain makes no sense: lower boundary exceeds upper boundary!')
    elseif size(dom,1)==1
        % Assume same domain for all vars
        dom = repmat(dom,nvars(1),1);
    elseif size(dom,1)~=nvars
        error('The number of specified domains should match the number of specified variables.')
    end
else
    error("The domain of the spatial variables should be specified as nx2 array of type 'double', for n variables.")
end

% % % Next, check which object in the PDE structure is actually declared
if strcmpi(type,'x') || strcmpi(type,'state')
    obj = 'x';
elseif strcmpi(type,'u') || strcmpi(type,'control') || strcmp(type,'actuate')
    obj = 'u';
elseif strcmpi(type,'w') || strcmpi(type,'input') || strcmpi(type,'in') || strcmpi(type,'disturbance')
    obj = 'w';
elseif strcmpi(type,'y') || strcmpi(type,'sense') || strcmpi(type,'observe') || strcmpi(type,'observation')
    obj = 'y';
elseif strcmpi(type,'z') || strcmpi(type,'output') || strcmp(type,'out') || strcmp(type,'regulate')
    obj = 'z';
else
    error('It is unclear what type of PDE object is desired; please specify one of ''state'', ''input'', ''control'', ''output'', or ''sense''.')
end


% % % Now then, declaring the actual object is not too exciting
% First, generate an index identifying this object;
%idx = getPDEvarID(obj);
idx = stateNameGenerator;

% Declare an object with the right size and variables.
out_var = pde_struct();
out_var.(obj){1}.size = sz;
out_var.(obj){1}.vars = vars;
out_var.(obj){1}.dom = dom;
out_var.(obj){1}.ID = idx;

% Set the object table.
out_var.([obj,'_tab']) = zeros(1,2);
out_var.([obj,'_tab'])(1,1) = idx;
out_var.([obj,'_tab'])(1,2) = sz;

% Set the variable as a free term to use for constructing an equation.      % DJ, 01/03/2025
out_var.free{1}.size = sz;
out_var.free{1}.vars = vars;
out_var.free{1}.term{1}.(obj) = 1;
out_var.free{1}.term{1}.C = eye(sz);
if strcmp(obj,'x')
    out_var.free{1}.term{1}.loc = vars';
    out_var.free{1}.term{1}.D = zeros(1,size(vars,1));
end
if strcmp(obj,'x') || strcmp(obj,'w') || strcmp(obj,'u')
    out_var.free{1}.term{1}.I = cell(nvars,1);
end

end


%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
function ID = getPDEvarID(obj)
% ID = GETPDEVARID(obj) generates a unique integer "ID" to identfy a PDE
% variables (state, input or output) of type "obj".
%
% INPUT
% - obj:    'char' belonging to either 'x', 'y', 'z', 'u', or 'w',
%           specifying what type of object the ID is generated for. This is
%           only so that e.g. a state and output variable can have the same
%           ID, as we can distinguish them from the fact that one is a
%           state and the other an output.
%
% OUTPUT
% - ID:     scalar integer of type 'double' specifying an index for a new
%           PDE variable of the specified type. The new ID corresponds to 1
%           plus the previous ID for that type of object. The counter
%           resets when Matlab is restarted.

persistent idx_x idx_u idx_w idx_y idx_z

    if strcmp(obj,'x')
        if isempty(idx_x)
            idx_x = 0;
        end
        idx_x = idx_x+1;
        ID = idx_x;
    elseif strcmp(obj,'u')
        if isempty(idx_u)
            idx_u = 0;
        end
        idx_u = idx_u+1;
        ID = idx_u;
    elseif strcmp(obj,'w')
        if isempty(idx_w)
            idx_w = 0;
        end
        idx_w = idx_w+1;
        ID = idx_w;
    elseif strcmp(obj,'y')
        if isempty(idx_y)
            idx_y = 0;
        end
        idx_y = idx_y+1;
        ID = idx_y;
    elseif strcmp(obj,'z')
        if isempty(idx_z)
            idx_z = 0;
        end
        idx_z = idx_z+1;
        ID = idx_z;
    end
end