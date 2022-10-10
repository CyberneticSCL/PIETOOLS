function out = convert(obj,convertTo)
arguments
    obj;
    convertTo {mustBeMember(convertTo,{'pie','ddf'})} = '';
end

if strcmp(convertTo, 'pie')
    out = sys();
    if strcmp(obj.type,'pde')
        out.params = convert_PIETOOLS_PDE(obj.params);
    elseif strcmp(obj.type,'dde')
        out.params = convert_PIETOOLS_DDE(obj.params,'pie');
    elseif strcmp(obj.type,'nds')
        out.params = convert_PIETOOLS_NDS(obj.params,'pie');
    elseif strcmp(obj.type,'ddf')
        out.params = convert_PIETOOLS_DDF(obj.params,'pie');
    end
    out.type = 'pie';
elseif strcmp(convertTo,'ddf') && (strcmp(obj.type,'nds')||strcmp(obj.type,'dde'))
    out = obj;
    if strcmp(obj.type,'nds')
        out.params = convert_PIETOOLS_NDS2DDF(obj.params);
    elseif strcmp(obj.type,'dde')
        out.params = minimize_PIETOOLS_DDE2DDF(obj.params);
    end
    out.type = 'ddf';
end

fprintf('Conversion to %s was successful\n', convertTo);
end