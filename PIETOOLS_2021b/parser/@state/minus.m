function sumTerms = minus(objA,objB)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that performs objA-objB
% Input: 
% objA, objB - state class objects
% Output:
% sumTerms - A*[objA;objB], where A=[1,-1]

if ~isa(objA,'state')||~isa(objB,'state')
    error('Only state type objects can be added together');
end

s.type = '.'; s.subs = 'veclength';

if sum(subsref(objA,s))~=sum(subsref(objB,s))
    error('States of unequal length cannot be added');
end

isdot_A = isdot(objA); isdot_B = isdot(objB); isout_A=isout(objA); isout_B = isout(objB);
s.type = '.'; s.subs = 'veclength';

if any((isdot_A&isout_B)|(isout_A&isdot_B)|(isout_A&isout_B)|(isdot_A&isdot_B))
    error("Any linear combinations of outputs and time-derivative of state are not allowed");
end

[objC,permMatsA,permMatsB] = combine(objA,objB); % objA = permMats{1}*objC and objB = permMats{2}*objC

% objC = state2terms(objC);

opvar T; 
T.R.R0 = permMatsA-permMatsB;

sumTerms = terms(T,objC);
end