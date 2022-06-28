function obj = vertcat(varargin)
if nargin==1
    obj = varargin{1};
else
    objA = varargin{1};
    objB = varargin{2};
    if isa(objA,'state')
        objA = state2terms(objA);
    end
    if isa(objB,'state')
        objB = state2terms(objB);
    end
    
    tempoperator = [objA.operator; objB.operator];
    tempstatevec = [objA.statevec; objB.statevec];
    
    obj = terms(tempoperator,tempstatevec);
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end