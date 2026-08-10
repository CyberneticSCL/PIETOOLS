function Cout = lrmultiply(L, C, R)
% Given C = A + B^T d, this function
% performs the operation 
% Cout = vec(L*C*R) 
% where C, Cout are structs with fields A, B^T
% and L/R are sparse matrices

A = C.A;
Bt = C.Bt;

% the following avoids construction of big (R'\otimes L)
X = reshape(A,size(L,2),size(R,1));
Y = L*X*R;  % (R^T\otimes L)*A = vec(L*X*R), X = reshape(A)
Aout = Y(:);

Btout = sparse(size(L,1)*size(R,2),size(Bt,2));
for i=1:size(Bt,2)
X = reshape(Bt(:,i),size(L,2),size(R,1));
Y = L*X*R;
Btout(:,i) = Y(:);
end

Cout = struct('A',Aout,'Bt', Btout);
end