function sumTerms = plus(objA,objB)
if ~isa(objA,'signals')||~isa(objB,'signals')
    error('Only state type objects can be added together');
end

s.type = '.'; s.subs = 'len';

if sum(subsref(objA,s))~=sum(subsref(objB,s))
    error('States of unequal length cannot be added');
end

[objC,permMatsA,permMatsB] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

% NOTE: need to replace the below with new opvar object
mi T; T.kernel = permMatsA+permMatsB;
sumTerms = termvar(T,objC);
end