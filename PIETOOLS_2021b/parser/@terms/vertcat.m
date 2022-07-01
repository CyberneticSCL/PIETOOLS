function obj = vertcat(varargin)
if nargin==1
    obj = varargin{1};
else
    objA = varargin{1};
    objB = varargin{2};
    if isempty(objA)
        obj = objB;
        return
    elseif isempty(objB)
        obj = objA;
        return
    end
    
    if isa(objA,'state')
        objA = state2terms(objA);
    end
    if isa(objB,'state')
        objB = state2terms(objB);
    end
    
    tempoperator = [objA.operator; objB.operator];
    tmpA = objA.statevec.state;tmpB = objB.statevec.state;
    tempstatevec.state = [tmpA;tmpB];
    tempstatevec.diff = [objA.statevec.diff; objB.statevec.diff];
    tempstatevec.delta = [objA.statevec.delta; objB.statevec.delta];
    
    obj = terms(tempoperator,tempstatevec);
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end