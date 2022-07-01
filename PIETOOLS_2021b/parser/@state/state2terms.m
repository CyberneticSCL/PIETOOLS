function out = state2terms(obj,operator,var,val)
arguments
    obj;
    operator {mustBeMember(operator,{'','diff','delta'})} = '';
    var = 0;
    val = 0;
end

opvar T; T.R.R0 = eye(length(obj));


termsOut = obj;

if strcmp(operator,'diff')
    diffvar =var;
    difforder = val;
    for i=1:length(termsOut)
        idx = find(isequal(termsOut(i).var,diffvar));
        termsOut(i).diff_order(idx) = difforder;
    end
elseif strcmp(operator,'delta')
    subvar = var;
    subval = val;
    for i=1:length(termsOut)
        idx = find(isequal(termsOut(i).var,subvar));
        termsOut(i).delta_val(idx) = subval;
    end
end

out = terms(T,termsOut);
end