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
    idx = find(isequal(obj(i).var,var));
    if ~strcmp("inf",obj(i).maxdiff)&&(obj(i).diffOrder(idx)+order>obj(i).maxdiff(idx))
        msg = "Differentiation with respect to "+var.varname+" exceeds specified max derivative, "+num2str(obj(i).maxdiff(idx));
        error(msg);
    else
        obj(i).diffOrder(idx) = obj(i).diffOrder(idx)+order;
    end
end
end