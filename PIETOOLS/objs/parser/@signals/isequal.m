function logval = isequal(A,B)
logval = 1;
if length(A)~=length(B)
    logval=0;
    return
end
if ~all(A.statename==B.statename)
    logval = 0;
    return
end
for i=1:length(A)
if ~all(A(i).diffOrder==B(i).diffOrder)
    logval=0;
    return
end
if ~all(A(i).maxdiff==B(i).maxdiff)
    logval=0;
    return
end
if ~all(A(i).dom==B(i).dom)
    logval=0;
    return
end
if ~isequal(A(i).var,B(i).var)
    logval = 0;
    return
end
if ~strcmp(A(i).type,B(i).type)
    logval=0;
    return
end
end

if ~all(A.len==B.len)
    logval = 0;
    return
end
if ~all(A.isInfinite==B.isInfinite)
    logval=0;
    return
end
end