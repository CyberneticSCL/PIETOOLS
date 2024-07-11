function out = int(obj, var, lim)
out = [];
for i=1:length(obj)
    tmpObj = obj(i);
    dom = tmpObj.dom;
    idx = find(isequal(tmpObj.var,var));
    if isempty(idx) % state drops out of integration
        K = int(1,var,lim(1),lim(2));
        out = [out;mtimes(K,tmpObj)];
    else
        K = eye([obj(i).len]);
        T = buildopvar('kernel',K,'lim',lim,'var',var,'dom',dom);
        out = [out;termvar(T,obj)];
    end
end
end