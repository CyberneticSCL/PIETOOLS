function [Aout,Bout] = lr_multiply(L,A,B,R)
% Given vec(C(d)) = A + B'*d, return the coefficients of vec(L*C(d)*R), using
% vec(L*C*R) = (R' kron L)*vec(C).
%
% The constant term is handled by reshaping, which avoids the Kronecker
% product entirely. For B the product is formed explicitly, as in
% @sdopvar/plus: here L and R are monomial selection matrices with a single
% nonzero per row, so the Kronecker product has only nnz(L)*nnz(R) nonzeros
% and applying it to all decision variables at once is far cheaper than
% looping over the rows of B, of which there may be many thousands.

X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;
Aout = Y(:);

Bout = B*kron(R.',L).';

end


