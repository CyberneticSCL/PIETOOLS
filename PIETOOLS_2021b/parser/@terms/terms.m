classdef terms
properties (SetAccess=protected)
    operator opvar;
    statevec struct;
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
    obj = plus(objA,objB);
    obj = minus(objA,objB);
    obj = uplus(objA);
    obj = uminus(objA);
    obj = diff(objA,var,order);
    obj = int(objA,var,lim);
    obj = delta(objA,var,val);
    obj = mtimes(objA,objB); % note only one of two inputs can be terms object
    obj = times(objA,objB); % note only one of two inputs can be terms object
    obj = horzcat(objA,objB);
    obj = vertcat(objA,objB);
end
end