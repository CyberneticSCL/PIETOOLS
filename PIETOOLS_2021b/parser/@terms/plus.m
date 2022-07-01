function sumTerms = plus(objA,objB)
if isa(objA,'state')
    objA = state2terms(objA);
elseif ~isa(objA,'terms')
    error('Only state/terms type objects can be added');
end
if isa(objB,'state')
    objB = state2terms(objB);
elseif ~isa(objB,'terms')
    error('Only state/terms type objects can be added');
end

if length(objA)~=length(objB)
    error('Terms of unequal length cannot be added');
end

tempoperator = [objA.operator objB.operator];
tmpA = objA.statevec.state;tmpB = objB.statevec.state;
tempstatevec.state = [tmpA;tmpB];
tempstatevec.diff = [objA.statevec.diff; objB.statevec.diff];
tempstatevec.delta = [objA.statevec.delta; objB.statevec.delta];
    
sumTerms = terms(tempoperator,tempstatevec);
end