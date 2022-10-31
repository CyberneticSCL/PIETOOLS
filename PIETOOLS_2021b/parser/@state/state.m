classdef (InferiorClasses={?polynomial,?dpvar})state
    properties
        type {mustBeMember(type,{'ode','pde','in','out'})} = 'ode';
        veclength {mustBeInteger,mustBePositive}=1;
        var {mustBeVector,mustBeA(var,["polynomial","double"])} = [pvar('t')];
        diff_order {mustBeInteger,mustBeVector,mustBeNonnegative}= [0];
    end
    properties (Hidden)
        maxdiff = "undefined";
        dom = [0,1];
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
                    obj.type = varargin{1};
                    if strcmp(varargin{1},'pde')
                        obj.var = [pvar('t'),pvar('s')];
                        obj.diff_order = [0,0];
                    end
                    obj.statename = stateNameGenerator();
                elseif nargin==2
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    if strcmp(obj.type,'pde')
                        obj.var = [pvar('t'),pvar('s')];
                        obj.diff_order = [0,0];
                    end
                    obj.statename = stateNameGenerator();
                elseif nargin==3
                    if size(varargin{3},1)~=1
                        error('var must be a row vector');
                    end
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    obj.var = varargin{3};
                    obj.statename = stateNameGenerator();
                    obj.diff_order = zeros(1,length(varargin{3}));
                elseif nargin==4 % internal use only, dont use this for constructing state vectors
                    obj.type = varargin{1};
                    obj.veclength = varargin{2};
                    obj.var = varargin{3};
                    obj.statename = varargin{4};
                    obj.diff_order = zeros(1,length(varargin{3}));
                elseif nargin>3
                    error('State class definition only takes 3 inputs');
                end
            end
        end
        
        function obj = set.dom(obj,dom)
            obj.dom = dom;
        end
        function obj = set.maxdiff(obj,maxdiff)
            obj.maxdiff = maxdiff;
        end

        % other class methods
        obj = subs(obj,var,var_val);
        obj = diff(obj,var,order);
        logval = eq(objA,objB);
        logval = isequal(objA,objB);
        obj = horzcat(varargin);
        obj = int(obj,var,limits);
        [logval,idx] = ismember(objA,objB)
        obj = minus(objA,objB);
        obj = mtimes(obj,K);
        logval = ne(objA,objB);
        obj = plus(objA,objB);
        obj = uplus(obj);
        obj = uminus(obj);
    end
end