function obj = removeequation(obj,eqnNumber)
obj.equation(eqnNumber) = [];
fprintf('%d equations were removed from the system\n', length(eqnNumber));
end