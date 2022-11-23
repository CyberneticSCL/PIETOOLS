function logval = isequal(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that tests objA=objB
% Input: 
% objA, objB - terms class objects 
% Output:
% logval - 0 or 1

logval = 1;
if ~isa(objA,'terms')||~isa(objB,'terms')
    logval = 0;
    return
end

if any(length(objA)~=length(objB))
    logval = 0;
    return
end

if (objA.operator~=objB.operator)
    logval=0;
    return
end
if ~isequal(objA.statevec,objB.statevec)
    logval=0;
    return
end
end