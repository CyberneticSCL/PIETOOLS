%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fcht.m     PIETOOLS 2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% If you modify this code, document all changes carefully and include date
% authorship, and a brief description of modifications
%
% Initial coding unknown

function a = fcht(v)
%FCHT	Fast Chebyshev transform
% FCHT(V) computes the Chebyshev transform of a N+1 by 1 array V.  If V
% corresponds to a function evaluated at the Chebyshev points
% cos(pi*(0:N)/N), then V is interpolated by a linear combinations of the
% Chebyshev polynomials with weights given by FCHT(V).
% 
% 
% Example:
% Approximate f(x) = exp(x) over [-1,1] as a linear combination of the
% first three Chebyshev polynomials.
% 
% x = cos(pi*(0:2)/2); % establish 3 Chebyshev grid points
% 
% V = exp(x); % evaluate f(x) at Chebyshev grid points
% 
% a = fcht(V);
% xx = linspace(-1,1); % create dense grid over domain
% g = a(1)*1 + a(2)*xx + a(3)*(2*xx.^2 - 1); % sum the first three Chebyshev
%               % polynomials with respect to their corresponding weights
% plot(xx,exp(xx),xx,g); % visualize the approximation
v = v(:);
N = length(v) - 1;
v = [v; flipud(v(2:N))];
a = real(fft(v))/N;
a = [a(1)/2; a(2:N); a(N+1)/2];
end