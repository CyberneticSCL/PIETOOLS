function obj = removeControl(obj,input)
loc = ismember([obj.states.statename],input.statename);
obj.ControlledInputs(find(loc))= 0;
tmpMsg = [num2str(length(input)), ' inputs were designated as disturbance inputs'];
disp(tmpMsg);
end