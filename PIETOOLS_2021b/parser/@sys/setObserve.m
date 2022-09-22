function obj = setObserve(obj,output)
loc = ismember(obj.states.statename,output.statename);
obj.ObservedOutputs(find(loc))= 1;
fprintf('%d outputs were designated as observed outputs\n', length(output));
end