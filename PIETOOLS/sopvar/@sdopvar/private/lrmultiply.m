function Cout = lrmultiply(L, C, R)
% Given C = A + B^T d, this function
% performs the operation 
% Cout = vec(L*C*R) 
% where C, Cout are structs with fields A, B
% and L/R are sparse matrices
%
% initial SS
% vectorization AT 09/08/26
%

A = C.A;
B = C.B;

% the following avoids construction of big (R'\otimes L)
X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;  % (R^T\otimes L)*A = vec(L*X*R), X = reshape(A)
Aout = Y(:);

Bout = sparse(size(B,1),size(L,1)*size(R,2));

% new implementation using kronecker product
% vec(L*X*R) = (R' otimes L)*vec(X) = vec(X)' (R otimes L') 
% Rkron = kron(speye(size(L,2)), R);
% Lkron = kron(L, speye(size(R,1)));
LRkron = kron(R, L');%;Lkron*Rkron;
Bout = B*LRkron;

% Previous implementation
% it was a loop over decision variables
% for i=1:size(B,1)  % row slicing is likely to me slow
%     X = reshape(B(i,:)',size(L,2),size(R,1));
%     Y = L*X*R;
%     Bout(i,:) = Y(:);
% end 
% if max(abs(Bout - Out), [], 'all') > 1.e-10
%     error('ERROR')
% end

Cout = struct('A',Aout,'B', Bout);
end