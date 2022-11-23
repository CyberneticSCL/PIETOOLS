function out = convert(obj,convertTo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sys-class method that converts obj from one type to another (typically
% PIE)
% Input: 
% obj - sys class object
% convertTo - string, 'pie' or 'ddf'
% Output:
% obj - sys class object with pie or ddf parameters

arguments
    obj;
    convertTo {mustBeMember(convertTo,{'pie','ddf'})} = '';
end

if strcmp(convertTo, 'pie')
    if strcmp(obj.type,'pde')
        tmp = convert_PIETOOLS_PDE(obj.params);
    elseif strcmp(obj.type,'dde')
        tmp = convert_PIETOOLS_DDE(obj.params,'pie');
    elseif strcmp(obj.type,'nds')
        tmp = convert_PIETOOLS_NDS(obj.params,'pie');
    elseif strcmp(obj.type,'ddf')
        tmp = convert_PIETOOLS_DDF(obj.params,'pie');
    elseif strcmp(obj.type,'pie')
        out = obj;
        return;
    end
    out = sys();
    out.type = 'pie';
    out.params = tmp;
elseif strcmp(convertTo,'ddf') && (strcmp(obj.type,'nds')||strcmp(obj.type,'dde'))
    out = obj;
    if strcmp(obj.type,'nds')
        tmp = convert_PIETOOLS_NDS2DDF(obj.params);
    elseif strcmp(obj.type,'dde')
        tmp = minimize_PIETOOLS_DDE2DDF(obj.params);
    end
    out.type = 'ddf';
    out.params = tmp;
end

fprintf('Conversion to %s was successful\n', convertTo);
end