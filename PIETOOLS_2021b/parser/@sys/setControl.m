function obj = setControl(obj,input)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that designates 'input' state object as a control input
% Input: 
% obj - sys class object
% input - state class object
% Output:
% obj - sys class object with 'input' designated as a control input

loc = ismember(obj.states.statename,input.statename);
obj.ControlledInputs(find(loc))= 1;
fprintf('%d inputs were designated as controlled inputs\n', length(input));
end