function obj = eq(objA,objB)
objA = +objA; % this ensures the objects are of termvar class
objB = +objB;
if ~isa(objA,'termvar')&&(objA==0) % if one of them is a zero double
    obj = -objB;
elseif ~isa(objB,'termvar')&&(objB==0)
    obj = objA;
elseif isa(objA,'termvar')&&isa(objB,'termvar')
    obj = objA-objB;
else
    error('Invalid equation format. Equations must be of the form "expr==0" or "exprA==exprB"');
end
end