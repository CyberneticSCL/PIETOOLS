function obj = diff(obj,var,order)
% Check inputs: if var is a string and if the order of differentiation is a valid integer
if ~isa(var,"polynomial")
    error('Variable of Differentiation must be a pvar object');
elseif ~isa(order,'double')||(order<0)||(floor(order)~=order)
    error('Order of Differentiation must be a non-negative integer');
end

% for a state object of the form [x; y; z;...] iterate through each state
for i=1:size(obj.len,2)
    tmpObj = obj(i);
    [loc,logval] = ismember(var,tmpObj.var); % find if var is present in obj(i)
    if logval
        tmp = zeros(size(tmpObj.diffOrder)); 
        tmp(loc) = order;
        tmpObj.diffOrder = tmpObj.diffOrder+tmp;
    else
        tmpObj.multiplier = diff(tmpObj.multiplier,var,order);
    end
    obj(i) = tmpObj;
end
end