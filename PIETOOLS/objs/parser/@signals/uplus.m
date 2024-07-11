function out = uplus(obj)
T = buildopvar('identity','dom',obj.dom);
out = termvar(T,obj);
end