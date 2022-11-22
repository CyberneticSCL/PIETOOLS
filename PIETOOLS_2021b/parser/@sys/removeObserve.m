function obj = removeObserve(obj,output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that designates 'output' state object as a regulated
% output
% Input: 
% obj - sys class object
% output - state class object
% Output:
% obj - sys class object with 'output' designated as a regulated output

loc = ismember(obj.states.statename,output.statename);
obj.ObservedOutputs(find(loc))= 0;
fprintf('%d outputs were designated as regulated outputs\n', length(output));
end