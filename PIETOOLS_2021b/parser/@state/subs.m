function obj = subs(obj,old,new)

if nargin~=3
    error('Substitution requires 3 inputs: state object, variable to be substituted, value to be substituted.');
end
if ~isa(new,'polynomial')&&((new~=0)&&(new~=1))
    error('Subs operation can only be performed with pvar variable or at the boundaries 0 and 1.');
end
for i=1:length(obj)
    if isequal(obj(i).var(1),old)
        error('Subs can only be performed on spatial variables');
    end
    if ismember(obj(i).type,{'ode','in','out'})
        error('Subs operation on ODE states, inputs, and outputs is not supported')
    end
    idx = find(isequal(obj(i).var,old));
    obj(i).var(idx) = new;
end
end