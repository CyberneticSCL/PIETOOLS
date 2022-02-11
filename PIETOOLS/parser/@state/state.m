classdef state 
    properties
        type {validateType(type)} = {'ode'};
        length {mustBeInteger,mustBePositive,mustBeScalarOrEmpty} = 1;
        var {validateVar(var)} = {[pvar('t')]};
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
                    if strcmp(varargin{1},'pde')
                        obj.var = {[pvar('t'),pvar('s1')]};
                    end
                elseif nargin==2
                    obj.type = varargin{1};
                    obj.length = varargin{2};
                    obj.segregation = 1;
                    if strcmp(obj.type,'pde')
                        obj.var = {[pvar('t'),pvar('s1')]};
                    end
                elseif nargin==3
                    if size(varargin{3},1)~=1
                        error('var must be a row vector');
                    end
                    obj.type = {varargin{1}};
                    obj.length = varargin{2};
                    obj.segregation = 1;
                    obj.var = {varargin{3}};
                elseif nargin==4 % internal use only, dont use this for constructing state vectors
                    obj.type = varargin{1};
                    obj.length = varargin{2};
                    obj.var = varargin{3};
                    obj.segregation = varargin{4};
                elseif nargin>4
                    error('State class definition only takes 3 inputs');
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
function validateVar(prop)
if ~iscell(prop)&&~isvector(prop)
    error('var must be column cell array of polynomial row vectors');
end
for i=1:length(prop)
    mustBeVector(prop{i});
    mustBeA(prop{i},'polynomial');
end
end