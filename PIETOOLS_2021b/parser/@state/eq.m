function terms_obj = eq(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that creates the equality constraint objA==objB
% Input: 
% objA, objB - state class objects of same length
% Output:
% terms_obj - an equation objA-objB==0

if ~isa(objA,'state')&&(objA==0)
    terms_obj = state2terms(-objB);
elseif ~isa(objB,'state')&&(objB==0)
    terms_obj = state2terms(objA);
elseif isa(objA,'state')&&isa(objB,'state')
    terms_obj = objA-objB;
else
    error('Invalid equation format. Equations must be of the form "expr==0" or "exprA==exprB"');
end
end