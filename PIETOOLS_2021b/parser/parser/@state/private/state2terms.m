function termsOut = state2terms(obj,operator,var,val)
arguments
    obj;
    operator {mustBeMember(operator,{'','diff','delta'})} = '';
    var = 0;
    val = 0;
end

termsOut = struct();
termsOut.state = obj;

termsOut.diff = cellfun(@(x) zeros(size(x)),obj.var,'un',0); 
termsOut.delta = obj.var;

if strcmp(operator,'diff')
    diffvar =var;
    difforder = val;
    idx = find(obj.var,diffvar);
    termsOut.diff(:,idx) = difforder;
elseif strcmp(operator,'delta')
    subvar = var;
    subval = val;
    idx = find(obj.var,subvar);
    termsOut.delta(:,idx) = subval;
end
end