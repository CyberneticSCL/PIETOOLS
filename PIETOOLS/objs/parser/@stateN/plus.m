function obj = plus(objA, objB)
if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be added together');
end
if length(objA)~=length(objB)
    error('States of unequal length cannot be added');
end

matA = eye(length(objA)); 
obj = vertcat(objA,objB); 
obj.multiplier = [matA, matA]*obj.multiplier;
obj = combine(obj);
end