function Pop = coeffs2opvar(C_str)
% POP = COEFFS2OPVAR(C_STR) takes a structure, C_STR, and builds an 'opvar'
% object with associated parameters, so that
%   Pop.R.R1 = Zd(s)'*C{1}*Zd(th)
% where s=C_str.var1 and th=C_str.var2

% Check the input
if ~isa(C_str,'struct')
    error("Input should be of type 'struct'.")
end
if ~isfield(C_str,'C')
    error("No coefficients specified.")
end
if ~isfield(C_str,'dim')
    error("Ambiguous dimensions of operator.")
end

% Initialize the operator
Pop = opvar();
if isfield(C_str,'dom')
    Pop.I = C_str.dom;
end
if isfield(C_str,'vars')
    Pop.var1 = C_str.vars(1);
    Pop.var2 = C_str.vars(2);
end
Pop.dim = [[0,0];C_str.dim];

m = C_str.dim(1);     n = C_str.dim(2);
C = C_str.C;
if ~isa(C,'cell') || ~numel(C)==3
    error("Coefficients should be specified as 1x3 cell");
end
C0 = C{1};  C1 = C{2};  C2 = C{3};
d = size(C1,1)/m - 1;
if any(size(C0)~=[m*(d+1),n]) || any(size(C1)~=size(C2)) || size(C1,2)~=n*(d+1)
    error("Coefficients should be specified as m*(d+1) x n*(d+1) arrays")
end

Z1 = polynomial(eye(d+1),(0:d)',Pop.var1.varname,[d+1,1]);
Z2 = polynomial(eye(d+1),(0:d)',Pop.var2.varname,[d+1,1]);

Pop.R.R0 = kron(eye(m),Z1')*C0;
Pop.R.R1 = kron(eye(m),Z1')*C1*kron(eye(n),Z2);
Pop.R.R2 = kron(eye(m),Z1')*C2*kron(eye(n),Z2);

end