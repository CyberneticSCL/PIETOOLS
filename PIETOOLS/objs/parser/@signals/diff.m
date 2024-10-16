function obj = diff(obj,var,order)
if nargin==1
    var = obj(1).var; var = var(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

if (length(var)>length(order))&&numel(order)==1
    order = order*ones(size(var));
end

for i=1:length(obj)
    for j=1:length(var)
        tmpvar = var(j);
        idx = find(isequal(obj(i).var,tmpvar));
        if ~strcmp("inf",obj(i).maxdiff(idx))&&(obj(i).diffOrder(idx)+order(j)>obj(i).maxdiff(idx))
            msg = "Differentiation with respect to "+tmpvar.varname+" exceeds specified max derivative, "+num2str(obj(i).maxdiff(idx));
            error(msg);
        else
            obj(i).maxdiff(idx) = obj(i).diffOrder(idx)+order(j);
        end
    end
end
end