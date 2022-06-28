function deltaTerms = delta(obj,var,val)
if nargin==1
    var = obj.var(1);
    val = 0;
elseif nargin==2
    val = 0;
elseif nargin>3
    error('Delta operator takes at most 3 inputs');
end

% start converting to terms object and then perform differentiation
deltaTerms = state2terms(obj,'delta',var,val);
end