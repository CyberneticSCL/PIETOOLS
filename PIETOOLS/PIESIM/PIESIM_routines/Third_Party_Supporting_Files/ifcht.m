%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ifcht.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding unknown

function v = ifcht(a)
%IFCHT	Fast Inverse Chebyshev transform
% IFCHT(V) computes the inverse Chebyshev transform of a N+1 by 1 array A.
% If A corresponds to the weights of a linear combination of Chebyshev
% polynomials, then IFCHT(a) computes the dataset interpolated by A at the
% Chebyshev points cos(pi*(0:N)/N).
% 
% 
% Example:
% Suppose A = [3; 2; 1].  Then the function 
%             f(x) = A(1)*1 + A(2)*x + A(3)*(2*x^2 - 1)
% evaluated at x = cos(pi*(0:2)/2) = [1,0,-1] is given by ifcht(A).
%
% x = cos(pi*(0:2)/2); %create Chebyshev grid of three points
% xx = linspace(-1,1); % create dense grid over domain
% 
% A = [3; 2; 1];
% f = A(1)*1 + A(2)*xx + A(3)*(2*xx.^2 - 1);
% 
% plot(xx,f,x,ifcht(A),'.','MarkerSize',20);
a = a(:);
N = length(a) - 1;
a = N*[a(1)*2; a(2:N); a(end)*2];
v = ifft([a; flipud(a(2:end-1))]);
v = v(1:N+1);
end