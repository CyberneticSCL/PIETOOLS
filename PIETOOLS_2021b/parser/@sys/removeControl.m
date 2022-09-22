function obj = removeControl(obj,input)
loc = ismember(obj.states.statename,input.statename);
obj.ControlledInputs(find(loc))= 0;
fprintf('%d inputs were designated as disturbance inputs\n', length(input));
end