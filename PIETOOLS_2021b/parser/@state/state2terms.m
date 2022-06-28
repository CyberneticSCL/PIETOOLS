function out = state2terms(obj,operator,var,val)
arguments
    obj;
    operator {mustBeMember(operator,{'','diff','delta'})} = '';
    var = 0;
    val = 0;
end

opvar T; T.R.R0 = eye(length(obj));

termsOut = struct();
termsOut.state = obj;

termsOut.diff = cellfun(@(x) zeros(size(x)),obj.var,'un',0); 
termsOut.delta = obj.var;

if strcmp(operator,'diff')
    diffvar =var;
    difforder = val;
    idx = cellfun(@(x) difforder*isequal(x,diffvar),obj.var,'un',0);
    termsOut.diff = idx;
elseif strcmp(operator,'delta')
    subvar = var;
    subval = val;
    idx = cellfun(@(x) subs(x,subvar,subval),obj.var,'un',0);
    termsOut.delta = idx;
end

out = terms(T,termsOut);
end