function logval = ne(objA,objB)
if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be compared');
end

logval = 1;

if all(objA.statename==objB.statename)
    logval=0;
end
end