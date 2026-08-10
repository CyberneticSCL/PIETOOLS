function Cout = lrmultiply(L, C, R)
% Given C = A + B^T d, this function
% performs the operation 
% Cout = vec(L*C*R) 
% where C, Cout are structs with fields A, B
% and L/R are sparse matrices

A = C.A;
B = C.B;

% the following avoids construction of big (R'\otimes L)
X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;  % (R^T\otimes L)*A = vec(L*X*R), X = reshape(A)
Aout = Y(:);

Bout = sparse(size(B,1),size(L,1)*size(R,2));
for i=1:size(B,1)  % row slicing is likely to me slow
X = reshape(B(i,:)',size(L,2),size(R,1));
Y = L*X*R;
Bout(i,:) = Y(:);
end

Cout = struct('A',Aout,'B', Bout);
end