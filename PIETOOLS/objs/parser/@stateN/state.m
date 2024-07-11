classdef (InferiorClasses={?polynomial,?dpvar})state
    properties (SetAccess=protected) 
        % default values on creation
        type = ["finite"; "ode"];
        len =[1];
        dom = {[]};
        maxdiff = {[inf]};
        var = {[pvar('t')]};
        multiplier = [1];
        intLim = {[]};
        diffOrder = {[0]};
    end
    properties (Hidden, SetAccess=protected)
        % for faster comparison of states; used internally
        statename;
    end
    methods (Access = {?sys,?state})
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
                %
                if nargin>7
                    obj.diffOrder = varargin{8};
                end
                if nargin>6
                    obj.intLim = varargin{7};
                end
                if nargin>5
                    obj.multiplier = varargin{6};
                end
                if nargin>4
                    obj.var = varargin{5};
                end
                if nargin>3
                    obj.maxdiff = varargin{4};
                end
                if nargin>2
                    obj.dom = varargin{3};
                end
                if nargin>1
                    obj.len = varargin{2};
                end
                if ~strcmp(varargin{1},'ode')
                    if length(varargin{1})==1
                        obj.type = {"infinite";varargin{1}};
                    else
                        obj.type = varargin{1};
                    end
                end
                if nargin==9 % internal use only; 
                    obj.statename = varargin{9};
                else
                    obj.statename = stateNameGenerator();
                end
            end
        end
        
        % THESE NEED TO BE MODIFIED FOR ROBUST INPUT CHECK
        function obj = set.type(obj,type)
            obj.type = type;
        end
        function obj = set.len(obj,len)
            obj.len = len;
        end
        function obj = set.var(obj,var)
            obj.var = var;
        end
        function obj = set.multiplier(obj,multiplier)
            obj.multiplier = multiplier;
        end
        function obj = set.intLim(obj,intLim)
            obj.intLim = intLim;
        end
        function obj = set.maxdiff(obj,maxdiff)
            obj.maxdiff = maxdiff;
        end
        function obj = set.dom(obj,dom)
            obj.dom = dom;
        end
        function obj = set.diffOrder(obj,diffOrder)
            obj.diffOrder = diffOrder;
        end

        % other class methods
        obj = diff(obj,var,order);
        obj = eq(objA,objB);
        obj = horzcat(varargin);
        obj = int(obj,var,limits);
        logval = isequal(objA,objB);
        [logval,idx] = ismember(objA,objB);
        val = length(obj);
        obj = minus(objA,objB);
        obj = mtimes(obj,K);
        logval = ne(objA,objB);
        obj = plus(objA,objB);
        obj = subs(obj,var,val);
        obj = subsref(obj,idx);
        obj = uminus(obj);
        obj = uplus(obj);
        obj = vertcat(varargin);
    end
end