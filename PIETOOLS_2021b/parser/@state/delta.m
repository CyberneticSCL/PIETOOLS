function deltaTerms = delta(obj,var,val)
if nargin==1
    var = obj.var(1);
    val = 0;
elseif nargin==2
    val = 0;
elseif nargin>3
    error('Delta operator takes at most 3 inputs');
end

if (val~=0)&& (val~=1)
    error('Dirac operation can only be performed at the boundaries 0 or 1');
end
if ismember(obj.type,{'ode','in','out'})
    error('Delta operation on ODE states, inputs, and outputs is not supported')
end

% start converting to terms object and then perform differentiation
deltaTerms = state2terms(obj,'delta',var,val);
end