function prodTerms = mtimes(K,objA)

if isa(objA,'state')&&isa(K,'state')
    error('Two state type objects cannot be multiplied');
end

if isa(K,'state')
    error('Left multiplication by a state object is not supported');
end

isdot_A = []; isout_A=[]; 
for i=1:length(objA)
    isdot_A = [isdot_A; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
    isout_A = [isout_A; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
end
if any((isdot_A|isout_A))
    error("Multiplication involving vectors with outputs or time-derivative of state is not allowed");
end


if numel(K)~=1 && (size(K,2)~=sum(objA.veclength))&& (sum(objA.veclength)~=1)
    error('Dimensions of multiplier and state vector do not match. Cannot be multiplied');
end

opvar T;
if numel(K)==1
    s.type = '.'; s.subs = 'veclength';
    T.R.R0 = K*eye(sum(subsref(objA,s)));
else
    T.R.R0 = K;
end
prodTerms = terms(T,objA);
end