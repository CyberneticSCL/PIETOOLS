% CLENCURT nodes x (Chebyshev points) and weights w
% for Clenshaw-Curtis quadrature
%
%   int_{-1}^{1} f(x) dx ~ sum_{k=0}^{N} w(k+1) f(x(k+1)),  x(k+1) = cos(pi*k/N).
%
% INPUT
% - N:  positive integer, the number of subintervals (N+1 nodes).
%
% OUTPUT
% - x:  (N+1) x 1 'double', the Chebyshev points cos(pi*(0:N)'/N), from 1
%       to -1;
% - w:  1 x (N+1) 'double', the quadrature weights.

% Written by Professor Lloyd Trefethen, Oxford University
% Obtained at http://people.maths.ox.ac.uk/trefethen/clencurt.m

function [x,w] = clencurt(N)
% Nodes x = cos(theta)
theta = pi*(0:N)'/N; x = cos(theta);
% Weights: end weights, then the interior weights 2*v/N from a cosine sum
% that differs for even and odd N
w = zeros(1,N+1); ii = 2:N; v = ones(N-1,1);
if mod(N,2)==0
w(1) = 1/(N^2-1); w(N+1) = w(1);
for k=1:N/2-1, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
v = v - cos(N*theta(ii))/(N^2-1);
else
w(1) = 1/N^2; w(N+1) = w(1);
for k=1:(N-1)/2, v = v - 2*cos(2*k*theta(ii))/(4*k^2-1); end
end
w(ii) = 2*v/N;
end
