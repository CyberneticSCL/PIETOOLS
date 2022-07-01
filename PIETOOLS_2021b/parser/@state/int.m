function prodTerms = int(objC, var, lim)

if any(isequal(lim,var))
    error('Limits of integration and integration variable are the same. Change the variable naming');
end

opvar T;
if poly2double(lim(2))
T.R.R2 = eye(sum(objC.veclength));
end
if poly2double(lim(1))
T.R.R1 = eye(sum(objC.veclength));    
end

for i=1:length(objC)
    idx = find(isequal(objC(i).var,var));
    if ~isempty(idx)
        objC(i).delta_val(idx) = pvar('theta');
    end
end

prodTerms = terms(T,objC);
end