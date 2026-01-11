function Pop = coeffs2opvar(C,dom,var1,var2)
% POP = COEFFS2OPVAR(C,VAR1,VAR2) takes a cell of coefficients, C1, and
% monomial degree, d, and builds an 'opvar' object with associated
% parameters, so that
%   Pop.R.R1 = Zd(s)'*C{1}*Zd(th)
% where s=var1 and th=var2


% Initialize the operator
Pop = opvar();
% Set the variables
if nargin>=2
    Pop.I = dom;
end
if nargin>=3
    if ispvar(var1)
        Pop.var1 = var1;
    elseif ischar(var1) || isstr(var1)
        Pop.var1 = pvar(var1);
    elseif iscellstr(var1)
        Pop.var1 = pvar(var1{1});
        if numel(var1)==2
            Pop.var2 = pvar(var1{2});
        elseif numel(var1)>2
            error("Too many variables specified")
        end
    else
        error("Variables should be specified as 'pvar' objects")
    end
end
if nargin>=4
    if ispvar(var2)
        Pop.var2 = var2;
    elseif ischar(var2) || isstr(var2)
        Pop.var2 = pvar(var2);
    elseif iscellstr(var2)
        Pop.var2 = pvar(var2{1});
        if numel(var2)>1
            error("Too many variables specified")
        end
    else
        error("Variables should be specified as 'pvar' objects")
    end
end

if ~isa(C,'cell') || ~numel(C)==2
    error("Coefficients should be specified as 1x2 cell");
end
C1 = C{1};  C2 = C{2};
if any(size(C1)~=size(C2)) || size(C1,1)~=size(C1,2)
    error("Coefficients should be specified as d+1 x d+1 arrays")
end
d = size(C1,1)-1;

Z1 = polynomial(eye(d+1),(0:d)',Pop.var1.varname,[d+1,1]);
Z2 = polynomial(eye(d+1),(0:d)',Pop.var2.varname,[d+1,1]);

Pop.R.R1 = Z1'*C1*Z2;
Pop.R.R2 = Z1'*C2*Z2;

end