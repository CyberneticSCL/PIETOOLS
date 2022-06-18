function prodTerms = mtimes(K,objA)

if isa(objA,'state')&&isa(K,'state')
    error('Two state type objects cannot be multiplied');
end

if isa(K,'state') && (size(objA,1)==1)
    prodTerms = mtimes(objA,K);
    return
end

objC = state2terms(objA);
opvar T; T.R.R0 = K;

prodTerms = terms(T,objC);
end