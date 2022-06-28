function prodTerms = mtimes(K,objA)

if isa(objA,'terms')&&isa(K,'terms')
    error('Two terms type objects cannot be multiplied');
end

if numel(K)~=1 && size(K,2)~=length(objA)
    error('Dimensions of multiplier and terms object do not match. Cannot be multiplied');
end

opvar T; 
if numel(K)==1
    T.R.R0 = K*eye(length(objA));
else
    T.R.R0 = K;
end

prodTerms = terms(T*objA.operator,objA.statevec);
end