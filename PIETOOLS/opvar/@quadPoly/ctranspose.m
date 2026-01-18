function H = ctranspose(F)
% This function performs conjugate transpose of a matrix
% valued polynomial.
%
% Input:
% F: a matrix-valued polynomial of size m-by-n in variables s, theta
%
% Output:
% H: conjugate transpose of F.

H = F;
m  = F.dim(1); n  = F.dim(2);
ds = size(F.Zs,1);
dt = size(F.Zt,1);

[I,J,V] = find(F.C);
if isempty(V)
    Cnew = sparse(n*ds, m*dt);
else
    a = mod(I-1,ds)+1;
    b = mod(J-1,dt)+1;
    i = floor((I-a)./ds)+1;
    j = floor((J-b)./dt)+1;

    Inew = (j-1)*ds + a;
    Jnew = (i-1)*dt + b;

    nRows = n*ds;

    if numel(V) > 2e5
        lin = Inew + (Jnew-1)*nRows;
        [~,p] = sort(lin);
        Inew = Inew(p); Jnew = Jnew(p); V = V(p);
    end

    Cnew = sparse(Inew, Jnew, conj(V), nRows, m*dt);
end

H = quadPoly(Cnew, F.Zs, F.Zt, [n m], F.ns, F.nt);
end
