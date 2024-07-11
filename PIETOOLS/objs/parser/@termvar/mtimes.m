function out = mtimes(K,objA)
if isa(K,'termvar')
    error('Left multiplication by state/terms object is not supported');
end
if numel(K)~=1 && size(K,2)~=length(objA)
    error('Dimensions of multiplier and terms object do not match. Cannot be multiplied');
end
if numel(K)==1
    K = K*eye(length(objA));
end
out = termvar(K*objA.operator,objA.state);
end