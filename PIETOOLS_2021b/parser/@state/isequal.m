function logval = isequal(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that specifies if objA=objB
% Input: 
% objA, objB - state class object
% Output:
% logval - 0 (if objA is not equal to objB) or 1 (if objA=objB)

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
if any(~isequal(objA(i).var,objB(i).var))
    logval=0;
    return
end
end
end