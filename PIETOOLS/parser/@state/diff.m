function terms_obj = diff(obj,var,order)
if nargin==1
    var = obj.var_indep(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

% start converting to terms object and then perform differentiation

end