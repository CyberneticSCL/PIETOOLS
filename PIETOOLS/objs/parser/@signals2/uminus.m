function out = uminus(obj)
T = buildopvar('identity','dom',obj.dom);
out = termvar(-T,obj);
end