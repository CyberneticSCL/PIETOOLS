function obj = removeequation(obj,eqnNumber)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that removes equation at location `eqnNumber' from list
% of equations
% Input: 
% obj - sys class object
% eqnNumber - positive integer scalar or vector
% Output:
% obj - sys class object with 'eqnNumber' removed

obj.equation(eqnNumber) = [];
fprintf('%d equations were removed from the system\n', length(eqnNumber));
end