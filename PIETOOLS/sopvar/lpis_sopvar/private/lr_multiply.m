function [Aout,Bout] = lr_multiply(L,A,B,R,K)
% function [Aout,Bout] = lr_multiply(L,A,B,R)                               % MMP, 09/26/2026 (was)
% Given vec(C(d)) = A + B'*d, return the coefficients of vec(L*C(d)*R), using
% vec(L*C*R) = (R' kron L)*vec(C).
%
% The constant term is handled by reshaping, which avoids the Kronecker
% product entirely. For B the product is formed explicitly, as in
% @sdopvar/plus: here L and R are monomial selection matrices with a single
% nonzero per row, so the Kronecker product has only nnz(L)*nnz(R) nonzeros
% and applying it to all decision variables at once is far cheaper than
% looping over the rows of B, of which there may be many thousands.
%
% Optional K = kron(R.',L).', for a caller applying the same L and R to     % MMP, 09/26/2026
% several coefficient sets (one per gamma cell): forming it once per pair   % MMP, 09/26/2026
% instead of once per cell removed ~3 of 11.7 s of a 2-D poscopvar build.   % MMP, 09/26/2026
%
% MMP, 09/26/2026: optional fifth input K, as above. Without it the
%                  behaviour is unchanged.

X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;
Aout = Y(:);

% Bout = B*kron(R.',L).';                                                   % MMP, 09/26/2026 (was)
if nargin<5 || isempty(K),  K = kron(R.',L).';  end                         % MMP, 09/26/2026
Bout = B*K;                                                                 % MMP, 09/26/2026

end
