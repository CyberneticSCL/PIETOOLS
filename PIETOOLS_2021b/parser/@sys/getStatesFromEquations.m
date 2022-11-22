function out = getStatesFromEquations(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that identifies unique states in equations list
% Input: 
% obj - sys class object
% Output:
% out - vector of state class objects

equations = obj.equation;
eqnNum = length(equations);
if eqnNum==0
    out = [];
else
out = [];
for i=1:eqnNum
    tempterms = equations(i);
    out = combine_statename(out,tempterms.statevec);
end
end
end
function out = combine_statename(varargin)
out = varargin{1};
for i=2:nargin
    tmp = varargin{i};
    for j=1:length(tmp)
        s.type = '()'; s.subs = {j};
        temp = subsref(tmp,s); 
        if isempty(out)
            out = temp;
        else
            if ~ismember(temp.statename,out.statename)
                out = [out; temp];
            end
        end
    end
end
end