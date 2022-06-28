function sumTerms = uplus(objA)
objC = state2terms(objA);
opvar T; 
T.R.R0 =eye(length(objA));
sumTerms = terms(T,objC.statevec);
end