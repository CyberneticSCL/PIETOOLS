function sumTerms = plus(objA,objB)
if ~isa(objA,'signals')||~isa(objB,'signals')
    error('Only state type objects can be added together');
end

if sum([objA.len])~=sum([objB.len])
    error('States of unequal length cannot be added');
end

[objC,permMatsA,permMatsB] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

% NOTE: need to replace the below with new opvar object
T = buildopvar('multiplier',permMatsA+permMatsB,'dom',objC.dom);
sumTerms = termvar(T,objC);
end