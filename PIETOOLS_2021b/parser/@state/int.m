function prodTerms = int(objA, var, lim)

if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end

objC = state2terms(objA);

opvar T; T.var2 = var; 
if poly2double(lim(2))
T.R.R2 = eye(length(objA));
end
if poly2double(lim(1))
T.R.R1 = eye(length(objA));    
end

prodTerms = terms(T,objC.statevec);
end