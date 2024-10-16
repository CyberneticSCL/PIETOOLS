function logval = isequal(objA,objB)
logval = 1;
if ~isa(objA,'termvar')||~isa(objB,'termvar')
    logval = 0;
    return
end
if ~isequal(objA.operator,objB.operator)
    logval=0;
    return
end
if ~isequal(objA.state,objB.state)
    logval=0;
    return
end
end