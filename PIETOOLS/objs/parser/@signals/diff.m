function obj = diff(obj,var,order)
if nargin==1
    tmpvar = obj(1).var;
    var = tmpvar{1}; var = var(1);
    order = 1;
elseif nargin==2
    order = 1;
elseif nargin>3
    error('Differentiation takes at most 3 inputs');
end

for i=1:length(obj)
    tmpvar = obj(i).var; tmpvar =tmpvar{1};
    idx = find(isequal(tmpvar,var));
    s1.type = '()'; s1.subs = {i};
    s2.type = '.'; s2.subs = 'diffOrder'; 
    s3.type = '.'; s3.subs = 'maxdiff';
    s4.type = '()'; s4.subs = {idx};
    smax = [s1;s3;s4]; sdiff = [s1;s2;s4];
    if ~strcmp("inf",subsref(obj, smax))&&(subsref(obj,sdiff)+order>subsref(obj,smax))
        msg = "Differentiation with respect to "+var.varname+" exceeds specified max derivative, "+num2str(subsref(obj,smax));
        error(msg);
    else
        obj(i).diffOrder(idx) = subsref(obj,sdiff)+order;
    end
end
end