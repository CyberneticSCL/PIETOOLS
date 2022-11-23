function obj = eq(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% terms-class method that creates the equality constraint objA==objB
% Input: 
% objA, objB - terms class objects
% Output:
% terms_obj - terms class object objA-objB

if isa(objA,'state')
    objA = state2terms(objA);
end
if isa(objB,'state')
    objB = state2terms(objB);
end
if ~isa(objA,'terms')&&(objA==0)
    obj = -objB;
elseif ~isa(objB,'terms')&&(objB==0)
    obj = objA;
elseif isa(objA,'terms')&&isa(objB,'terms')
    obj = objA-objB;
else
    error('Invalid equation format. Equations must be of the form "expr==0" or "exprA==exprB"');
end
end