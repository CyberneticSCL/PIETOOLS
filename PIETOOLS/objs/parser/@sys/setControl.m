function obj = setControl(obj,input)
loc = ismember([obj.states.statename],input.statename);
obj.ControlledInputs(find(loc))= 1;
tmpMsg = [num2str(length(input)), ' inputs were designated as controlled inputs'];
disp(tmpMsg);
end