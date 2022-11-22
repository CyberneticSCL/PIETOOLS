function out = getParams(obj)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that parses equations into standard parametric form
% Input: 
% obj - sys class object
% Output:
% out - sys class object with parameters initialized

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
