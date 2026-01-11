function Pop = coeffs2opvar2d(C_str)
% POP = COEFFS2OPVAR2d(C_STR) takes a structure, C_STR, and builds an 
% 'opvar2d' object with associated parameters, so that
%   Pop.R22{i,j} = Zd(s)'*C_str.C{i,j}*Zd(th)
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
Pop = opvar2d();
if isfield(C_str,'dom')
    Pop.I = C_str.dom;
end
if isfield(C_str,'vars')
    Pop.var1 = C_str.vars(:,1);
    Pop.var2 = C_str.vars(:,2);
end
Pop.dim = [zeros(3,2);C_str.dim];

m = C_str.dim(1);     n = C_str.dim(2);
C = C_str.C;
if ~isa(C,'cell') || ~numel(C)==9
    error("Coefficients should be specified as 3x3 cell");
end
%C0 = C{1};  C1 = C{2};  C2 = C{3};
d = C_str.deg;
% d = size(C1,1)/m - 1;
% if any(size(C0)~=[m*(d+1),n]) || any(size(C1)~=size(C2)) || size(C1,2)~=n*(d+1)
%     error("Coefficients should be specified as m*(d+1) x n*(d+1) arrays")
% end

Z11 = polynomial(eye(d+1),(0:d)',Pop.var1(1).varname,[d+1,1]);
Z12 = polynomial(eye(d+1),(0:d)',Pop.var1(2).varname,[d+1,1]);
Z21 = polynomial(eye(d+1),(0:d)',Pop.var2(1).varname,[d+1,1]);
Z22 = polynomial(eye(d+1),(0:d)',Pop.var2(2).varname,[d+1,1]);

Z1 = kron(Z11,Z12);
Z2 = kron(Z21,Z22);

Pop.R22{1,1} = kron(eye(m),Z1')*C{1,1};
Pop.R22{2,1} = kron(eye(m),Z1')*C{2,1}*kron(eye(n),Z21);
Pop.R22{3,1} = kron(eye(m),Z1')*C{3,1}*kron(eye(n),Z21);
Pop.R22{1,2} = kron(eye(m),Z1')*C{1,2}*kron(eye(n),Z22);
Pop.R22{1,3} = kron(eye(m),Z1')*C{1,3}*kron(eye(n),Z22);
Pop.R22{2,2} = kron(eye(m),Z1')*C{2,2}*kron(eye(n),Z2);
Pop.R22{3,2} = kron(eye(m),Z1')*C{3,2}*kron(eye(n),Z2);
Pop.R22{2,3} = kron(eye(m),Z1')*C{2,3}*kron(eye(n),Z2);
Pop.R22{3,3} = kron(eye(m),Z1')*C{3,3}*kron(eye(n),Z2);

end