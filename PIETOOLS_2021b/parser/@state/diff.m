function obj = diff(obj,var,order)
if nargin==1
    var = obj.var(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

for i=1:length(obj)
    if isequal(var,obj(i).var(1))&&((order>1)||(any(obj(i).diff_order>0)))
        error('Higher order derivatives of time and temporal-spatial mixed derivatives are currently not supported')
    end
    if ismember(obj(i).type,{'in','out'})
        error('Differentiation of inputs and outputs is currently not supported')
    end
    idx = find(isequal(obj(i).var,var));
    obj(i).diff_order(idx) = obj(i).diff_order(idx)+order;
end
end