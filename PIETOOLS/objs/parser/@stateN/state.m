classdef (InferiorClasses={?polynomial,?dpvar})state
    properties (SetAccess=protected)
        type = ["finite"; "ode"];
        len =[1];
        var = {[pvar('t')]};
        multiplier = [1];
        intLim = {[]};
        diffOrder = {[0]};
        maxdiff = {[inf]};
        dom = {[]};
    end
    properties (Hidden, SetAccess=protected)
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
            end
        end
        
        function obj = set.type(obj,type)
        end
        function obj = set.len(obj,len)
        end
        function obj = set.var(obj,var)
        end
        function obj = set.multiplier(obj,multiplier)
        end
        function obj = set.intLim(obj,intLim)
        end
        function obj = set.maxdiff(obj,maxdiff)
        end
        function obj = set.dom(obj,dom)
        end
        function obj = set.diffOrder(obj,diffOrder)
        end

        % other class methods
        obj = diff(obj,var,order);
        obj = eq(objA,objB);
        obj = horzcat(varargin);
        obj = vertcat(varargin);
        obj = int(obj,var,limits);
        obj = minus(objA,objB);
        obj = mtimes(obj,K);
        obj = plus(objA,objB);
        obj = subs(obj,var,val);
        obj = subsref(obj,idx);
        obj = uminus(obj);
        obj = uplus(obj);
        logval = isequal(objA,objB);
        [logval,idx] = ismember(objA,objB);
        logval = ne(objA,objB);
    end
end