function [x,w] = gaussleg(n,a,b)
% Gauss-Legendre nodes and weights on [a,b], exact for degree <= 2n-1.
% Golub-Welsch on the Jacobi matrix: stable at the orders used here.
k = 1:n-1;
beta = k./sqrt(4*k.^2-1);
[V,Dg] = eig(diag(beta,1)+diag(beta,-1));
[x,ix] = sort(diag(Dg));
w = 2*(V(1,ix).^2)';
x = 0.5*(b-a)*x + 0.5*(a+b);
w = 0.5*(b-a)*w;
end
