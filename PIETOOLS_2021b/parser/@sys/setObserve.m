function obj = setObserve(obj,output)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that designates 'output' state object as an observed
% output
% Input: 
% obj - sys class object
% output - state class object
% Output:
% obj - sys class object with 'output' designated as an observed output

loc = ismember(obj.states.statename,output.statename);
obj.ObservedOutputs(find(loc))= 1;
fprintf('%d outputs were designated as observed outputs\n', length(output));
end