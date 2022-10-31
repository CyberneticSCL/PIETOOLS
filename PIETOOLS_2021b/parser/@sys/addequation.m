function obj = addequation(obj,eqn)
if isa(eqn,'state')
    eqn = state2terms(eqn);
elseif isa(eqn,'terms')
    obj.equation = [obj.equation; eqn];
    fprintf('%d equations were added to sys() object\n',length(eqn));
else
    error('Unknown equation type. Cannot be added to system.');
end
end