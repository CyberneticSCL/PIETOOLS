function obj = addequation(obj,eqn)
if isa(eqn,'signals')
    eqn = +eqn;
end
if isa(eqn,'termvar')
    obj.equation = [obj.equation; eqn];
    tmpMsg = [num2str(length(eqn))+" equations were added to sys() object"];
    disp(tmpMsg);
else
    tmpMsg = 'Unknown equation type. Cannot be added to system.'; 
    error(tmpMsg);
end
end