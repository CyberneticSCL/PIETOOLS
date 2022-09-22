function obj = removeObserve(obj,output)
loc = ismember(obj.states.statename,output.statename);
obj.ObservedOutputs(find(loc))= 0;
fprintf('%d outputs were designated as regulated outputs\n', length(output));
end