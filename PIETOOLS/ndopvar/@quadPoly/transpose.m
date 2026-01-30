function H = transpose(F)
% TRANSPOSE Transpose of quadPoly WITHOUT swapping left/right variables.
%
% If F(s,t) = (I_m ⊗ Zs(s)^T) * C * (I_n ⊗ Zt(t)),
% then H(s,t) = F(s,t)' is represented as
%   H(s,t) = (I_n ⊗ Zs(s)^T) * Cnew * (I_m ⊗ Zt(t)),
% using the SAME Zs/ns on the left and SAME Zt/nt on the right.
%
% Cnew is obtained by transposing only the outer (m×n) block structure of C,
% while keeping the within-block (ds×dt) coordinates fixed.

m  = F.dim(1);
n  = F.dim(2);

ds = prod(cellfun(@numel, F.Zs));  % tensor basis size on s-side
dt = prod(cellfun(@numel, F.Zt));  % tensor basis size on t-side

C = sparse(F.C);
[i, j, v] = find(C);

% Within-block indices and block indices
a  = mod(i-1, ds) + 1;        % 1..ds
ib = floor((i-1) / ds);       % 0..m-1

b  = mod(j-1, dt) + 1;        % 1..dt
jb = floor((j-1) / dt);       % 0..n-1

% Swap ONLY the outer block indices (ib <-> jb)
i2 = jb * ds + a;             % rows: n blocks of size ds
j2 = ib * dt + b;             % cols: m blocks of size dt

Cnew = sparse(i2, j2, v, n*ds, m*dt);

H = quadPoly(Cnew, F.Zs, F.Zt, [n, m], F.ns, F.nt);
end