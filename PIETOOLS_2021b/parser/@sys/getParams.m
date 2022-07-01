function out = getParams(obj)
if strcmp(obj.type,'pde')
    obj.params = getPDEparams(obj);
    out = obj;
elseif strcmp(obj.type,'dde')
    error('Symbolic definition of DDEs is not currently supported and will be added in the future');
elseif strcmp(obj.type,'ddf')
    error('Symbolic definition of DDFs is not currently supported and will be added in the future');
elseif strcmp(obj.type,'nds')
    error('Symbolic definition of NDSs is not currently supported and will be added in the future');
end
end
