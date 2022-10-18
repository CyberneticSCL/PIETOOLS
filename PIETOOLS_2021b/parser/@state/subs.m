function obj = subs(obj,old,new)

if nargin~=3
    error('Substitution requires 3 inputs: state object, variable to be substituted, value to be substituted.');
end
if ~isa(new,'polynomial')&&((new~=0)&&(new~=1))
    error('Subs operation can only be performed with pvar variable or at the boundaries 0 and 1.');
end
for i=1:length(obj)
    if isequal(obj(i).var(1),old)
        if ~poly2double(obj(i).var(1)-new)
            error('Subs on time can only be performed to add delays: must be of the form t-tau, where tau is positive real');
        end
        if poly2double(obj(i).var(1)-new) && double(obj(i).var(1)-new)<0
            error('Positive delays are not allowed since the system is non-causal');
        end
    end
%     if ismember(obj(i).type,{'ode','in','out'})
%         error('Subs operation on ODE states, inputs, and outputs is not supported')
%     end
    idx = find(isequal(obj(i).var,old));
    obj(i).var(idx) = new;
end
end