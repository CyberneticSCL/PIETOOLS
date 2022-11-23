function logval = isdot(objA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that specifies if objA has time derivatives in it.
% Input: 
% objA - state class object
% Output:
% logval - logical array specifying if vector, objA, has components with time
% derivative

s.type = '.'; s.subs = 'veclength';
logval = []; 
for i=1:length(objA)
    logval = [logval; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
end
logval = boolean(logval);
end