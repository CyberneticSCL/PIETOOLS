function out = getStatesFromEquations(obj)
equations = obj.equation;
eqnNum = length(equations);
if eqnNum==0
    out = [];
else
out = equations{1}.statevec.state;
for i=2:eqnNum
    out = combine(out,equations{i}.statevec.state);
end
end