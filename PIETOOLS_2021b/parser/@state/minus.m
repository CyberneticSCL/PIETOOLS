function sumTerms = minus(objA,objB)

if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be added together');
end

if length(objA)~=length(objB)
    error('States of unequal length cannot be added');
end

for i=1:length(objA)
    if (objA(i).diff(1))&&(objB(i).diff(1))
        error("Two states with time derivatives cannot be added together");
    end
end

[objC,permMatsA,permMatsB] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

% objC = state2terms(objC);

opvar T; 
T.R.R0 = permMatsA-permMatsB;

sumTerms = terms(T,objC);
end