function obj = removeObserve(obj,output)
loc = ismember(obj.states.statename,output.statename);
obj.ObservedOutputs(find(loc))= 0;
end