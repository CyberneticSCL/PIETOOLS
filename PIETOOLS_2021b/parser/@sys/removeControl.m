function obj = removeControl(obj,input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that designates 'input' state object as a disturbance
% input
% Input: 
% obj - sys class object
% input - state class object
% Output:
% obj - sys class object with 'input' designated as a disturbance

loc = ismember(obj.states.statename,input.statename);
obj.ControlledInputs(find(loc))= 0;
fprintf('%d inputs were designated as disturbance inputs\n', length(input));
end