function prodTerms = mtimes(K,objA)
if isa(objA,'signals')&&isa(K,'signals')
    error('Two state type objects cannot be multiplied');
end
if isa(K,'signals')
    error('Left multiplication by a state object is not supported');
end

s.type = '.'; s.subs = 'len';
if numel(K)~=1 && (size(K,2)~=sum(subsref(objA,s)))&& (sum(subsref(objA,s))~=1)
    error('Dimensions of multiplier and state vector do not match. Cannot be multiplied');
end

mi T;
if numel(K)==1
    s.type = '.'; s.subs = 'veclength';
    T.kernel = K*eye(sum(subsref(objA,s)));
else
    T.kernel = K;
end
prodTerms = termvar(T,objA);
end