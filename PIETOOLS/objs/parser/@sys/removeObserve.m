function obj = removeObserve(obj,output)
loc = ismember([obj.states.statename],output.statename);
obj.ObservedOutputs(find(loc))= 0;
tmpMsg = [num2str(length(output)), ' outputs were designated as regulated outputs'];
disp(tmpMsg);
end