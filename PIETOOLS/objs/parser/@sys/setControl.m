function obj = setControl(obj,input)
inputID = [input.statename];
stateID = [obj.states.statename];
if ~all(strcmp(input.type,'in'))
    tmpMsg = ['Some or all arguments passed are not signals of type "input". Cannot be designated as controlled inputs.'];
    error(tmpMsg);
end
inputsAdded = 0;
for i=1:length(inputID)
    if ismember(inputID,stateID)&&~ismember(inputID, obj.CInames)
        obj.CInames = [obj.CInames; inputID(i)];
        inputsAdded = inputsAdded + input(i).len;
    end
end
tmpMsg = [num2str(inputsAdded), ' inputs were designated as controlled inputs'];
disp(tmpMsg);
end