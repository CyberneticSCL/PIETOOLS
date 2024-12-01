function out = int(obj, var, lim)
intOp = buildopvar('dom',{obj.dom},'var',{obj.var},'lim',lim, 'intvar',var);
out = termvar(intOp,obj);
end