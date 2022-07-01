function obj = setControl(obj,input)
loc = ismember(obj.states.statename,input.statename);
obj.ControlledInputs(find(loc))= 1;
end