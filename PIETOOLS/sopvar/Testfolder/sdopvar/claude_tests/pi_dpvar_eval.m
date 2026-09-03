function P = pi_dpvar_eval(D,dvars,dval)
% P = PI_DPVAR_EVAL(D,DVARS,DVAL) evaluates a 'dpvar' object at the given
% values of its decision variables, returning the resulting 'polynomial'.
%
% The 'dpvar' class provides 'subs' for independent variables only, so this
% fills in the missing operation. It is used to compare operators returned by
% 'poslpivar' against those returned by 'possopvar'.
%
% INPUTS
% - D:      'dpvar' object;
% - dvars:  cell array of decision variable names indexing DVAL;
% - dval:   array of values, with DVAL(i) the value of DVARS{i};
%
% OUTPUTS
% - P:      'polynomial' object of the same dimensions as D;
%
% MP, 08/22/2026: Initial coding

if ~isa(D,'dpvar')
    error("Input must be a 'dpvar' object.")
end
dvars = cellstr(string(dvars(:)));
dval = reshape(dval,[],1);

[m,n] = size(D);
dnames = cellstr(string(D.dvarname(:)));
nd = numel(dnames);
nZ = size(D.degmat,1);

% Weight of each row within a block: the first row holds the constant term,
% the remaining rows the coefficients of each decision variable.
if nd==0
    w = 1;
else
    [tf,loc] = ismember(dnames,dvars);
    if ~all(tf)
        error("Decision variable '"+string(dnames{find(~tf,1)})+"' has no assigned value.")
    end
    w = [1;dval(loc)];
end

% Contract the decision variable rows of each block, leaving one polynomial
% coefficient per monomial and matrix entry. The 'polynomial' constructor
% expects the coefficients of entry (p,q) in column p+m*(q-1).
coef = zeros(nZ,m*n);
for p = 1:m
    rows = (p-1)*(nd+1)+(1:nd+1);
    Cp = D.C(rows,:);
    for q = 1:n
        coef(:,p+m*(q-1)) = (w.'*Cp(:,(q-1)*nZ+(1:nZ))).';
    end
end

if isempty(D.varname)
    P = polynomial(reshape(sum(coef,1),m,n));
else
    P = polynomial(coef,D.degmat,D.varname,[m,n]);
end

end
