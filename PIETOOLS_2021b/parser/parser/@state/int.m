function prodTerms = int(objA, var, lim)

objC = state2terms(objA);

opvar T; T.var2 = var;
if ~isa(lim(2),'polynomial')
T.R.R2 = eye(objA.length);
end
if ~isa(lim(1),'polynomial')
T.R.R1 = eye(objA.length);    
end

prodTerms = terms(T,objC);
end