classdef (InferiorClasses={?state})terms
properties (SetAccess=protected)
    operator opvar;
    statevec state;
end
methods
    function obj = terms(varargin)
        if nargin>2
            error('Terms object constructor takes only two inputs');
        end
        if nargin>1
            obj.statevec = varargin{2};
        end
        if nargin>0
            obj.operator = varargin{1};    
        end
    end
    
%     obj = delta(objA,var,val);
%     obj = diff(objA,var,order);
%     obj = horzcat(varargin);
    obj = int(objA,var,lim);
    obj = length(obj);
    obj = minus(objA,objB);
    obj = eq(objA,objB);
    obj = isequal(objA,objB);
    obj = mtimes(objA,objB); % note only one of two inputs can be terms object
    obj = plus(objA,objB);
%     obj = times(objA,objB); % note only one of two inputs can be terms object
    obj = uplus(objA);
    obj = uminus(objA);
    obj = vertcat(varargin);
end
end