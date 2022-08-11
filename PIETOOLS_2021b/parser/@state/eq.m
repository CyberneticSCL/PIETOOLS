function logval = eq(objA,objB)
logval = 1;
if ~isa(objA,'state')||~isa(objB,'state')
    logval = 0;
    return
end

if any(size(objA)~=size(objB))
    logval = 0;
    return
end

for i=1:size(objA,1)
if (objA(i).statename~=objB(i).statename)
    logval=0;
    return
end
if any(objA(i).diff_order~=objB(i).diff_order)
    logval=0;
    return
end
end
end