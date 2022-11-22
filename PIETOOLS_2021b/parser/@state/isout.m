function logval = isout(objA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that specifies if objA has outputs in it.
% Input: 
% objA - state class object
% Output:
% logval - logical array specifying if vector, objA, has components with
% outputs

s.type ='.'; s.subs = 'veclength';
logval=[];
for i=1:length(objA)
    logval = [logval; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
end
logval = boolean(logval);
end