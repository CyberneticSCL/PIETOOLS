function obj = removeObserve(obj,output)
obj.OOnames(ismember(obj.OOnames,[output.statename]))= 0;
tmpMsg = [num2str(sum([output.len])), ' outputs were designated as regulated outputs'];
disp(tmpMsg);
end