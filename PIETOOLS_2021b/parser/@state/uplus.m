function sumTerms = uplus(objA)
s.type = '.'; s.subs = 'veclength';
opvar T; T.R.R0 = eye(sum(subsref(objA,s)));
sumTerms = terms(T,objA);
end