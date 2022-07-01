function diffTerms = diff(obj,var,order)
if nargin==1
    var = obj.var(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

if (val>1)&& (var==obj.var(1))
    error('Second and higher order derivatives of time is currently not supported')
end
if ismember(obj.type,{'in','out'})
    error('Differentiation of inputs and outputs is currently not supported')
end

% start converting to terms object and then perform differentiation
diffTerms = state2terms(obj,'diff',var,order);
end