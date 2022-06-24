function diffTerms = diff(obj,var,order)
if nargin==1
    var = obj.var(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

% start converting to terms object and then perform differentiation
objC = state2terms(obj,'diff',var,order);

opvar T; T.R.R0 = eye(length(obj));
diffTerms = terms(T,objC);
end