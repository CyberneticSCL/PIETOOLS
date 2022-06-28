function prodTerms = int(objA, var, lim)

if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end

try
    intterms = double(objA.operator.R.R2)|double(objA.operator.R.R1);
    if any(intterms(:))
        error('Double integrals are currently not supported')
    end
catch
    error('Double integrals are currently not supported')
end


opvar T; T.var2 = var;
if poly2double(lim(2))
T.R.R2 = objA.operator.R.R0;
end
if poly2double(lim(1))
T.R.R1 = objA.operator.R.R0;    
end

prodTerms = terms(T,objA.statevec);
end