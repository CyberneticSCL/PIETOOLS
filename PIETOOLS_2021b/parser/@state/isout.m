function logval = isout(objA)
logval=[];
for i=1:length(objA)
    logval = [logval; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
end
logval = boolean(logval);
end