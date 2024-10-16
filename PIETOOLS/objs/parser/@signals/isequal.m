function logval = isequal(A,B)
logval = 1;
if length(A)~=length(B)
    logval=0;
    return
end
for i=1:length(A)
if ~all(A(i).statename==B(i).statename)
    logval = 0;
    return
end
if ~all(isequal(A(i).diffOrder,B(i).diffOrder))
    logval=0;
    return
end
if ~all(isequal(A(i).maxdiff,B(i).maxdiff))
    logval=0;
    return
end
if ~all(isequal(A(i).dom,B(i).dom))
    logval=0;
    return
end
tmpA = A(i).var; tmpB = B(i).var;
if ~isequal(tmpA{1},tmpB{1}) % pvar objects have a weird problem
    logval = 0;
    return
end
if ~strcmp(A(i).type,B(i).type)
    logval=0;
    return
end
if ~all(A(i).len==B(i).len)
    logval = 0;
    return
end
end
end