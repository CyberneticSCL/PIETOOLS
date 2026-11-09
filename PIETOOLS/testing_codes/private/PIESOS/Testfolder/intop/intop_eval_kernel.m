function K = intop_eval_kernel(U,k,vals,d0)
% K = INTOP_EVAL_KERNEL(U,K,VALS,D0) evaluates the kth kernel of the
% 'intop' object U at the dummy variable values VALS and decision variable
% values D0, directly from the definition of the class,
%   K(s) = kron(I_m,[1,d']) * params{k} * kron(I_n,z(s)),
% with z(s) the Kronecker product of the monomials defined by Z.in.
%
% INPUTS
% - U:      m x n 'intop' object;
% - k:      index of the parameter to evaluate;
% - vals:   containers.Map from the names of the dummy variables
%           ('s1_dum1',...) to numeric values;
% - d0:     nd x 1 array of values for the decision variables U.Zd;
%
% OUTPUTS
% - K:      m x n array, the kernel evaluated at the specified values;
%
% NOTES
% - Independent of the 'kernels' method and of the constructor from
%   kernels, so usable as ground truth for both.

% DJ (CC), 10/15/2026: Initial coding;

[m,n] = size(U);
nd = numel(U.Zd);
vn = U.varnames.in_dum;
N = numel(vn);
% Monomials: Kronecker product, first variable slowest
z = 1;
for i=1:N
    Zi = U.Z.in{i};
    zi = ones(size(Zi,1),1);
    for kk=1:size(Zi,2)
        sval = vals(vn{i}{kk});
        zi = zi.*sval.^Zi(:,kk);
    end
    z = kron(z,zi);
end
P = U.params{k};
sel = kron(speye(m),[1,reshape(d0,1,[])]);
Cd = full(sel*P);                       % m x n*nZ
nZ = numel(z);
K = zeros(m,n);
for j=1:n
    K(:,j) = Cd(:,(j-1)*nZ+1:j*nZ)*z;
end
if nd~=numel(d0)
    error("Number of decision variable values must match number of decision variables.")
end

end
