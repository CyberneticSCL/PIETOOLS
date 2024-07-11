function obj = vertcat(varargin)
if nargin==1
    obj = varargin{1};
else
    objA = varargin{1};
    objB = varargin{2};
    
    objA = +objA; objB = +objB;
    if ~isa(objA,'termvar') || ~isa(objB,'termvar')
        error('Only state/terms type objects can be added');
    end
    
    if isempty(objA)
        obj = objB;
        return
    elseif isempty(objB)
        obj = objA;
        return
    end
    
   tempoperator = blkdiag(objA.operator,objB.operator);
   tempstate = vertcat(objA.state, objB.state);
    
    obj = terms(tempoperator,tempstate);
    if nargin>2
        obj = vertcat(obj,varargin{3:end});
    end
end
end