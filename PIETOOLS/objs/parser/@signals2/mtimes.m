function prodTerms = mtimes(K,objA)
if isa(K,'signals')
    error('Left multiplication by a state object is not supported');
end

if numel(K)~=1 && (size(K,2)~=sum([objA.len]))&& (sum([objA.len])~=1)
    error('Dimensions of multiplier and state vector do not match. Cannot be multiplied');
end

if numel(K)==1
    K = K*eye(sum([objA.len]));
end
T = buildopvar('multiplier',K,'dom',objA.dom);
prodTerms = termvar(T,objA);
end