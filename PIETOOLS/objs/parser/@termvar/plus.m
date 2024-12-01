function out = plus(objA,objB)
objA = +objA; objB = +objB;
if ~isa(objA,'termvar') || ~isa(objB,'termvar')
    error('Only state/terms type objects can be added');
end

if length(objA)~=length(objB)
    error('Terms of unequal length cannot be added');
end

tempoperator = [objA.operator objB.operator];
tempstate = vertcat(objA.state, objB.state);

out = termvar(tempoperator,tempstate);
end