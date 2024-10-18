function obj = removeequation(obj,eqnNumber)
obj.equation(eqnNumber) = [];
tmpMsg = [num2str(length(eqnNumber)), ' equations were removed from the system'];
disp(tmpMsg);
end