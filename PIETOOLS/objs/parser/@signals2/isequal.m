function logval = isequal(A,B)
logval = 1;
if length(A)~=length(B)
    logval=0;
    return
end
if ~all(isequal([A.statename],[B.statename]))
    logval = 0;
    return
end
if ~all(isequal([A.diffOrder],[B.diffOrder]))
    logval=0;
    return
end
if ~all(isequal([A.len],[B.len]))
    logval = 0;
    return
end

for i=1:length(A)
if ~all(isequal(A(i).maxdiff,B(i).maxdiff))
    logval=0;
    return
end
if ~all(isequal(A(i).dom,B(i).dom))
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
end