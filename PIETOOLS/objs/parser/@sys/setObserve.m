function obj = setObserve(obj,output)
outputID = [output.statename];
stateID = [obj.states.statename];
if ~all(strcmp(output.type,'out'))
    tmpMsg = ['Some or all arguments passed are not signals of type "output". Cannot be designated as observed outputs.'];
    error(tmpMsg);
end
outputsAdded = 0;
for i=1:length(outputID)
    if ismember(outputID,stateID)&&~ismember(outputID, obj.CInames)
        obj.OOnames = [obj.OOnames; outputID(i)];
        outputsAdded = outputsAdded + output(i).len;
    end
end
tmpMsg = [num2str(outputsAdded), ' outputs were designated as observed outputs'];
disp(tmpMsg);
end