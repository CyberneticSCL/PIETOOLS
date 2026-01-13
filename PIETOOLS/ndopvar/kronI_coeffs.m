function D_str = kronI_coeffs(C_str,q)
% D_STR = KRONI_COEFFS(C_STR,Q) takes the kronecker product of the operator 
% defined by C_STR with an identity matrix of dimension Q
%
% INPUTS
% - C_str:  struct representing a PI operator Pop;
% - q:      scalar integer specifying the dimension of the identity
%           operator;
%
% OUTPUTS
% - D_str:  struct representing the operator (Pop o I_{q});


% Construct the permutation matrices (I_{m(d+1)} o P)^T and 
% (I_{n(d+1)} o P) for P defined such that
%   Z_{d}(s) o I_{q} = P*(I_{q} o Z_{d}(s))
d = C_str.deg;
m = C_str.dim(1);   n = C_str.dim(2);

r_idcs = 1:q*(d+1);
c_idcs = (0:q-1)'*(d+1) + (1:d+1);
P = sparse(r_idcs,c_idcs,1,q*(d+1),q*(d+1));
Pt = sparse(c_idcs,r_idcs,1,q*(d+1),q*(d+1));
IP_m = spIkron(m,Pt);
IP_n = spIkron(n,P);

% Construct the coefficients representing C_str o I_{q}
C_cell = C_str.C;
D_cell = cell(size(C_cell));

D_cell{1} = IP_m*kron(C_cell{1},speye(q));
for ii=2:3
    D_cell{ii} = IP_m*kron(C_cell{ii},speye(q))*IP_n;
end
D_str = C_str;
D_str.dim = q*C_str.dim;
D_str.C = D_cell;


end