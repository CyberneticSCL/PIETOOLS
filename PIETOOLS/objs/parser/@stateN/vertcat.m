function obj = vertcat(varargin)
objA = varargin{1};
objB = varargin{2};

obj.type = [objA.type, objB.type];
obj.len = [objA.len, objB.len];
obj.var = [objA.var, objB.var];
obj.multiplier = [objA.multiplier, zeros(length(objA),length(objB)); zeros(length(objB),length(objA)), objB.multiplier];
obj.intLim = [objA.intLim, objB.intLim];        
obj.diffOrder = [objA.diffOrder, objB.diffOrder];
obj.maxdiff = [objA.maxdiff, objB.maxdiff];
obj.dom = [objA.dom, objB.dom];

obj.statename = [objA.statename, objB.statename];

if nargin>2
    obj = vertcat(obj,varargin{3:end});
end
end