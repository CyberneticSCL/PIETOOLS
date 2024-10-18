function obj = setObserve(obj,output)
loc = ismember([obj.states.statename],output.statename);
obj.ObservedOutputs(find(loc))= 1;
tmpMsg = [num2str(length(output)), ' outputs were designated as observed outputs'];
disp(tmpMsg);
end