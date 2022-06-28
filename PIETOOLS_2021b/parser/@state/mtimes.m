function prodTerms = mtimes(K,objA)

if isa(objA,'state')&&isa(K,'state')
    error('Two state type objects cannot be multiplied');
end

if isa(K,'state') && (numel(objA)==1)
    prodTerms = mtimes(objA,K);
    return
elseif isa(K,'state')
    error('Multiplication with a state object should have the state object on the right side of the operator');
end

if numel(K)~=1 && size(K,2)~=length(objA)
    error('Dimensions of multiplier and state vector do not match. Cannot be multiplied');
end

objC = state2terms(objA);
opvar T;
if numel(K)==1
    T.R.R0 = K*eye(length(objA));
else
    T.R.R0 = K;
end
prodTerms = terms(T,objC.statevec);
end