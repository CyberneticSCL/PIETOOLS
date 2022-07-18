function prodTerms = mtimes(K,objA)

if isa(objA,'state')&&isa(K,'state')
    error('Two state type objects cannot be multiplied');
end

if isa(K,'state')
    error('Left multiplication by a state object is not supported');
end

if numel(K)~=1 && (size(K,2)~=sum(objA.veclength))&& (sum(objA.veclength)~=1)
    error('Dimensions of multiplier and state vector do not match. Cannot be multiplied');
end

opvar T;
if numel(K)==1
    s.type = '.'; s.subs = 'veclength';
    T.R.R0 = K*eye(subsref(objA,s));
else
    T.R.R0 = K;
end
prodTerms = terms(T,objA);
end