function obj = addequation(obj,eqn)
if isa(eqn,'state')
    eqn = state2terms(eqn);
end

obj.equation = [obj.equation; eqn];
end