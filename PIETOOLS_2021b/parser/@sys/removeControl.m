function obj = removeControl(obj,input)
loc = ismember(obj.states.statename,input.statename);
obj.ControlledInputs(find(loc))= 0;
end