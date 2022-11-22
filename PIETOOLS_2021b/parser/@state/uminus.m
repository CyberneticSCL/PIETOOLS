function sumTerms = uminus(objA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that converts objA to terms format
% Input: 
% objA - state class objects
% Output:
% out - A*[objA], where A is an negative-identity PI operator

isdot_A = isdot(objA); isout_A=isout(objA); 
if any((isdot_A|isout_A))
    error("Unitary minus involving vectors with outputs or time-derivative of state is not allowed");
end
s.type = '.'; s.subs = 'veclength';
opvar T; T.R.R0 = -eye(sum(subsref(objA,s)));
sumTerms = terms(T,objA);
end