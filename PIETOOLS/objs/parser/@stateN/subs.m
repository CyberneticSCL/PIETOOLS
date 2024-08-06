function obj = subs(obj,old,new)
for i=1:length(obj.var)
    tmp = obj(i).var;
    idx= find(isequal(tmp,old));
    tmp(idx) = new;
    obj(i).var = tmp;        
end
end