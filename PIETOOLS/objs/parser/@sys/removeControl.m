function obj = removeControl(obj,input)
obj.CInames(ismember(obj.CInames,[input.statename]))= 0;
tmpMsg = [num2str(sum([input.len])), ' inputs were designated as disturbance inputs'];
disp(tmpMsg);
end