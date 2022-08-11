function prodTerms = mtimes(K,objA)

if isa(K,'terms')
    error('Left multiplication by state/terms object is not supported');
end
if isa(objA,'terms')&&isa(K,'terms')
    error('Two terms type objects cannot be multiplied');
end

isdot_A = isdot(objA); isout_A=isout(objA); 
% isdot_A = []; isout_A=[]; 
% for i=1:length(objA)
%     isdot_A = [isdot_A; objA(i).diff_order(1)*ones(subsref(objA(i),s),1)];
%     isout_A = [isout_A; strcmp(objA(i).type,'out')*ones(subsref(objA(i),s),1)];
% end
if any((isdot_A|isout_A))
    error("Multiplication involving vectors with outputs or time-derivative of state is not allowed");
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