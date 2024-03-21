function obj = mtimes(K,objA)
if isa(K,'state')
    error('Left multiplication by a state object is not supported');
end
obj = objA;
if length(K)==1
    K = K*eye(length(objA));
end
if length(objA)~=size(K,2)
    error('Inner dimensions of the multiplier and state do not match');
end
obj.multiplier = K*objA.multiplier;
obj = combine(obj);
end