function logval = isequal(objA,objB)
logval = 1;
if ~isa(objA,'state')||~isa(objB,'state')
    logval = 0;
    return
end

if length(objA)~=length(objB)
    logval = 0;
    return
end

for i=1:size(objA.len,1)
if (objA(i).statename~=objB(i).statename)
    logval=0;
    return
end
if any(objA(i).type~=objB(i).type)
    logval=0;
    return
end
if any(~isequal(objA(i).var,objB(i).var))
    logval=0;
    return
end
if any(objA(i).multiplier~=objB(i).multiplier)
    logval=0;
    return
end
if any(objA(i).intLim~=objB(i).intLim)
    logval=0;
    return
end
if any(objA(i).diffOrder~=objB(i).diffOrder)
    logval=0;
    return
end
if any(objA(i).maxDiff~=objB(i).maxDiff)
    logval=0;
    return
end
if any(objA(i).dom~=objB(i).dom)
    logval=0;
    return
end
end
end