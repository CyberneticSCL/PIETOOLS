function prodTerms = int(objA, var, lim)

objC = state2terms(objA);

opvar T; T.var2 = var;
if ~isa(lim(2),'polynomial')
T.R.R2 = eye(length(objA));
end
if ~isa(lim(1),'polynomial')
T.R.R1 = eye(length(objA));    
end

prodTerms = terms(T,objC);
end