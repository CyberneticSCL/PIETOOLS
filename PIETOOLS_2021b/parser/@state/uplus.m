function sumTerms = uplus(objA)
isdot_A = isdot(objA); isout_A=isout(objA); 
if any((isdot_A|isout_A))
    error("Unitary plus involving vectors with outputs or time-derivative of state is not allowed");
end

s.type = '.'; s.subs = 'veclength';
opvar T; T.R.R0 = eye(sum(subsref(objA,s)));
sumTerms = terms(T,objA);
end