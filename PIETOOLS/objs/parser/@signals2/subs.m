function obj = subs(obj,old,new)
if nargin~=3
    error('Substitution requires 3 inputs: state object, variable to be substituted, value to be substituted.');
end
if (length(old)>length(new))&&length(new)==1
    new = new*ones(size(old));
end

for i=1:length(obj)
    for j=1:length(old)
        idx = find(isequal(obj(i).var,old(j)));
        obj(i).var(idx) = new(j);
    end
end
end