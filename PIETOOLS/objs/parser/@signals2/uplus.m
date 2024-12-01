function out = uplus(obj)
T = buildopvar('multiplier',eye(sum([obj.len])),'dom',obj.dom);
out = termvar(T,obj);
end