function obj = subs(obj,old,new)
if nargin~=3
    error('Substitution requires 3 inputs: state object, variable to be substituted, value to be substituted.');
end
if ~isa(new,'polynomial')||~isa(new,"double")
    error('Subs operation can only be performed with pvar variable or double.');
end
for i=1:length(obj)
    idx = find(isequal(obj(i).var,old));
    obj(i).var(idx) = new;
end
end