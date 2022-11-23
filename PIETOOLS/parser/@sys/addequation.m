function obj = addequation(obj,eqn)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that adds eqn to list of equations in obj
% Input: 
% obj - sys class object
% eqn - terms class object
% Output:
% obj - sys class object with an new equation in the list

if isa(eqn,'state')
    eqn = state2terms(eqn);
elseif isa(eqn,'terms')
    obj.equation = [obj.equation; eqn];
    fprintf('%d equations were added to sys() object\n',length(eqn));
else
    error('Unknown equation type. Cannot be added to system.');
end
end