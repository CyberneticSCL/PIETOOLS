function out = uminus(objA)
isdot_A = isdot(objA.statevec); isout_A=isout(objA.statevec); 
if any((isdot_A|isout_A))
    error("Unitary minus involving vectors with outputs or time-derivative of state is not allowed");
end
objA.operator = -objA.operator;
out = objA;
end