classdef state 
    properties
        type {validateType(type)} = {'ode'};
        vecLength {mustBeInteger,mustBePositive,mustBeScalarOrEmpty} = 1;
        diff_order {validateDiff(diff_order)} = {[1]};
        var_indep {validateVar(var_indep)} = {[pvar('t')]};
    end
    properties (Hidden, SetAccess=protected)
        segregation {mustBeInteger,mustBePositive,mustBeVector}=[1];
    end
    methods
        function obj = state(varargin) %constructor
            if nargout==0
                for i=1:nargin
                    assignin('caller', varargin{i}, state());
                end
            else
                if nargin==1
                    obj.type = {varargin{1}};
                    if strcmp(varargin{1},'ode')
                        obj.diff_order = {[1]};
                    elseif strcmp(varargin{1},'pde')
                        obj.diff_order = {[1,0]};
                        obj.var_indep = {[pvar('t'),pvar('s1')]};
                    else
                        obj.diff_order = {[0]};
                    end
                elseif nargin==2
                    obj.type = varargin{1};
                    obj.vecLength = varargin{2};
                    obj.segregation = 1;
                    if strcmp(obj.type,'ode')
                        obj.diff_order = {[1]};
                    elseif strcmp(obj.type,'pde')
                        obj.diff_order = {[1,0]};
                        obj.var_indep = {[pvar('t'),pvar('s1')]};
                    else
                        obj.diff_order = {[0]};
                    end
                elseif nargin==3
                    if size(varargin{3},1)~=1
                        error('diff_order must be a row vector');
                    end
                    obj.type = {varargin{1}};
                    obj.vecLength = varargin{2};
                    obj.segregation = 1;
                    if strcmp(varargin{1},'pde')
                        if numel(varargin{3})==1
                            obj.diff_order = {[1,varargin{3}]};
                        else
                            obj.diff_order = {varargin{3}};
                        end
                        n_vars = numel(varargin{3});
                        vars(1) = pvar('t');
                        for i=2:n_vars
                            vars(i) = pvar(['s',num2str(i-1)]);
                        end
                        obj.var_indep = {vars};
                    elseif numel(varargin{3})==1
                       obj.diff_order = {varargin{3}};
                    else
                       error('Non-PDE State objects are functions of time and cannot have array input for order of differentiability');
                    end
                elseif nargin==4
                    if size(varargin{3},1)~=1
                        error('diff_order must be a row vector');
                    end
                    if size(varargin{4},1)~=1
                        error('var_indep must be a row vector');
                    end
                    obj.type = {varargin{1}};
                    obj.vecLength = varargin{2};
                    obj.segregation = 1;
                    if numel(varargin{3})==1
                        obj.diff_order = {[1,varargin{3}]};
                    else
                        obj.diff_order = {varargin{3}};
                    end
                    obj.var_indep = {varargin{4}};
                    obj = fix_property_dim(obj);
                elseif nargin==5 % internal use only, dont use this for constructing state vectors
                    obj.type = varargin{1};
                    obj.vecLength = varargin{2};
                    obj.diff_order = varargin{3};
                    obj.var_indep = varargin{4};
                    obj.segregation = varargin{5};
                elseif nargin>5
                    error('State class definition only takes 4 inputs');
                end
            end
        end
        
        % other class methods
        obj = diff(obj,var,order);
        obj = delta(obj,var,var_val);
        obj = plus(obj1,obj2);
        obj = minus(obj1,obj2);
        obj = uplus(obj);
        obj = uminus(obj);
        obj = mtimes(obj,K);
        obj = int(obj,var,limits);
        obj = horzcat(varargin);
        obj = vertcat(varargin);
    end
end
function validateType(prop)
if ~iscell(prop)&&~isvector(prop)
    error('Type must be cell column array');
end
if ~all(ismember(prop,{'ode','pde','out','in'}))
    error("Type must be one of the following strings: 'ode','pde','out','in'");
end
end
function validateDiff(prop)
if ~iscell(prop)&&~isvector(prop)
    error('diff_order must be cell column array');
end
for i=1:length(prop)
    mustBeVector(prop{i});
    mustBeInteger(prop{i});
    mustBeNonnegative(prop{i});
end
end
function validateVar(prop)
if ~iscell(prop)&&~isvector(prop)
    error('var_indep must be cell column array');
end
for i=1:length(prop)
    mustBeVector(prop{i});
    mustBeA(prop{i},'polynomial');
end
end