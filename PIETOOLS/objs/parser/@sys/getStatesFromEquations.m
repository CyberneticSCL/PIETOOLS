function out = getStatesFromEquations(obj)
equations = obj.equation;
eqnNum = length(equations);
if eqnNum==0
    out = [];
else
out = [];
for i=1:eqnNum
    tempterms = equations(i);
    out = combine_statename(out,tempterms.state);
end
end
end
function out = combine_statename(varargin)
out = varargin{1};
for i=2:nargin
    tmp = varargin{i};
    for j=1:length(tmp)
        temp = tmp(j); 
        if isempty(out)
            out = temp;
        else
            if ~ismember([temp.statename],out.statename)
                out = [out; temp];
            end
        end
    end
end
end