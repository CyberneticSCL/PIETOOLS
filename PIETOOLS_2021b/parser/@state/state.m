classdef state
    properties
        type {validateType(type)} = {'ode'};
        veclength {mustBeInteger,mustBePositive,mustBeVector}=[1];
        var {validateVar(var)} = {[pvar('t')]};
    end
    properties (Hidden, SetAccess=protected)
        statename;
    end
    methods (Access = {?terms, ?sys, ?state})
        objterms = state2terms(obj,operator,var,val);
        [out, varargout] = combine(varargin)
    end
    
    methods
        function obj = state(varargin) %constructor
            if nargout==0
                for i=1:nargin
                    obj = state();
                    obj.statename = stateNameGenerator();
                    assignin('caller', varargin{i}, obj);
                end
            else
                if nargin==1
                    obj.type = {varargin{1}};
                    if strcmp(varargin{1},'pde')
                        obj.var = {[pvar('t'),pvar('s1')]};
                    end
                    obj.statename = stateNameGenerator();
                elseif nargin==2
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    if strcmp(obj.type,'pde')
                        obj.var = {[pvar('t'),pvar('s1')]};
                    end
                    obj.statename = stateNameGenerator();
                elseif nargin==3
                    if size(varargin{3},1)~=1
                        error('var must be a row vector');
                    end
                    obj.type = {varargin{1}};
                    obj.veclength = varargin{2};
                    obj.var = {varargin{3}};
                    obj.statename = stateNameGenerator();
                elseif nargin==4 % internal use only, dont use this for constructing state vectors
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    obj.var = varargin{3};
                    obj.statename = varargin{4};
                elseif nargin>3
                    error('State class definition only takes 3 inputs');
                end
            end
        end
        
        % other class methods
        obj = delta(obj,var,var_val);
        obj = diff(obj,var,order);
        obj = eq(obj1,obj2);
        obj = horzcat(varargin);
        obj = int(obj,var,limits);
        obj = length(obj);
        obj = minus(obj1,obj2);
        obj = mtimes(obj,K);
        obj = ne(obj1,obj2);
        obj = plus(obj1,obj2);
        obj = size(obj);
        obj = uplus(obj);
        obj = uminus(obj);
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