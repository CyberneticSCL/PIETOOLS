function obj = diff(obj,var,order)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State-class method that differentiates a vector of states, x, in variable
% 'var' by 'order' times
% Input: 
% obj - state class object, x 
% var - pvar w.r.t. which differentiation is done
% order - positive integer, specifying order of differentiation
% Output:
% obj - state class object, \partial_{var}^{order} x

if nargin==1
    var = obj.var(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

for i=1:length(obj)
%     if isequal(var,obj(i).var(1))&&((order>1)||(any(obj(i).diff_order>0)))
%         error('Higher order derivatives of time and temporal-spatial mixed derivatives are currently not supported')
%     end
    if ismember(obj(i).type,{'in','out'})
        error('Differentiation of inputs and outputs is currently not supported')
    end
    idx = find(isequal(obj(i).var,var));
    if ~strcmp("undefined",obj(i).maxdiff)&&(obj(i).diff_order(idx)+order>obj(i).maxdiff(idx))
        msg = "Differentiation with respect to "+var.varname+" exceeds specified max derivative, "+num2str(obj(i).maxdiff(idx));
        error(msg);
    else
        obj(i).diff_order(idx) = obj(i).diff_order(idx)+order;
    end
end
end