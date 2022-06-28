function sumTerms = plus(objA,objB)
if isa(objA,'state')
    objA = state2terms(objA);
elseif ~isa(objA,'terms')
    error('Only state/terms type objects can be added')
end
if isa(objB,'state')
    objB = state2terms(objB);
elseif ~isa(objB,'terms')
    error('Only state/terms type objects can be added')
end

if length(objA)~=length(objB)
    error('Terms of unequal length cannot be added');
end

sumTerms = terms([objA.operator objB.operator],[objA.statevec;objB.statevec]);
end