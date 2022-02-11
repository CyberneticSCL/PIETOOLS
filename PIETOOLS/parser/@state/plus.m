function sumTerms = plus(objA,objB)

if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be added together');
end

opvar T1 T2;

[objC,permMats] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

T1.P = eye(finiteVecLength); T1.R.R0 = eye(infiniteVeclLength);
T2.P = eye(finiteVecLength2); T2.R.R0 = eye(infiniteVeclLength2);
T1 =T1*permMats{1}; T2 = T2*permMats{2};

objC = state2terms(objC);
sumTerms = terms(T1+T2,objC);
end