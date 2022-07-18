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
    
    if isempty(objA)
        obj = objB;
        return
    elseif isempty(objB)
        obj = objA;
        return
    end
    
    opvar zeroAB zeroBA;
    zeroAB.dim = [0 0; objA.operator.dim(2,1) objB.operator.dim(2,2)]; 
    zeroBA.dim = [0 0; objB.operator.dim(2,1) objA.operator.dim(2,2)];
    
    tempoperator = [objA.operator zeroAB; zeroBA objB.operator];
    tempstatevec = vertcat(objA.statevec, objB.statevec);
    
    obj = terms(tempoperator,tempstatevec);
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end