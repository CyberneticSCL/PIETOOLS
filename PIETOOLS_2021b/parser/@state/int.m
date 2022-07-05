function prodTerms = int(objC, var, lim)

if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end

if ~poly2double(lim(2))&&~poly2double(lim(1))
    error('Integral should have at least one limit of integration as 0 or 1');
end

opvar T;
if poly2double(lim(2))
    if (double(lim(2))==1)
        T.R.R2 = eye(sum(objC.veclength));
    else 
        error('Upper limit of integration can only be a "pvar" variable or 1');
    end
end
if poly2double(lim(1))  
    if (double(lim(1))==0)
        T.R.R1 = eye(sum(objC.veclength));    
    else
        error('Lower limit of integration can only be a "pvar" variable or 0');
    end
end
for i=1:length(objC)
    idx = find(isequal(objC(i).var,var));
    if ~isempty(idx)
        objC(i).delta_val(idx) = pvar('theta');
    end
end
prodTerms = terms(T,objC);
end