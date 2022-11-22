function logval = ne(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that specifies if objA=objB
% Input: 
% objA, objB - state class object
% Output:
% logval - 0 (if objA is equal to objB) or 1 otherwise

logval = ~(isequal(objA,objB));
end