function mu = inh_fdmax(a1,a2,b1,b2,c,N)
% MU = INH_FDMAX(A1,A2,B1,B2,C,N) largest real part of the eigenvalues of
% the second-order central finite-difference discretization, on the N x N
% interior grid of [0,1]^2 (h = 1/(N+1), homogeneous Dirichlet), of
% L x = a1 x_s1s1 + a2 x_s2s2 + b1 x_s1 + b2 x_s2 + c x; coefficients are
% function handles of (s1,s2). For N with a2 = [] the 1-D operator on [0,1]
% (a1, b1, c handles of s1). SCRATCH, 11/16/2026.
h = 1/(N+1);    s = (1:N)'*h;   e = ones(N,1);  I = speye(N);
D2 = spdiags([e -2*e e],-1:1,N,N)/h^2;  D1 = spdiags([-e 0*e e],-1:1,N,N)/(2*h);
if isempty(a2)
    L = spdiags(a1(s),0,N,N)*D2 + spdiags(b1(s),0,N,N)*D1 + spdiags(c(s),0,N,N);
else
    [S1,S2] = ndgrid(s,s);  dg = @(f) spdiags(f(S1(:),S2(:)),0,N^2,N^2);
    L = dg(a1)*kron(I,D2) + dg(a2)*kron(D2,I) + dg(b1)*kron(I,D1) + dg(b2)*kron(D1,I) + dg(c);
end
ev = eigs(L,4,'largestreal');   mu = max(real(ev));
end
