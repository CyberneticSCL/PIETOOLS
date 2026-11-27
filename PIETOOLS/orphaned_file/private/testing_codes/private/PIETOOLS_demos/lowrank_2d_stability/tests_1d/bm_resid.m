function [F,J,q] = bm_resid(w,P,rv)
% Preconditioned residual  F = W*(A(YY') - b)  and its analytic Jacobian.
% With the preconditioner of bm_setup, ||F|| is the Frobenius distance from
% Q = blkdiag(Y_i Y_i') to the affine set {A(X) = b}.
% w = [ z ; Y_1(:) ; ... ; Y_B(:) ],  Y_i is Ns(i) x rv(i)  (rv(i)=0 -> Q_i=0).
B = numel(P.Ns);
q = zeros(P.Ntot,1);
k = P.Kf;
if P.Kf>0, q(1:P.Kf) = w(1:P.Kf); end
Y = cell(1,B);
for i=1:B
    N = P.Ns(i); r = rv(i);
    if r==0, Y{i} = zeros(N,0); continue; end
    Y{i} = reshape(w(k+(1:N*r)),N,r);  k = k+N*r;
    Qi = Y{i}*Y{i}';
    q(P.rows{i}) = Qi(:);
end
F = P.W*(P.Ssym*q - P.bf);
if nargout<2, return; end
Jc = cell(1,B+1);
if P.Kf>0, Jc{1} = P.W*P.Ssym(:,1:P.Kf); else, Jc{1} = zeros(P.mres,0); end
for i=1:B
    N = P.Ns(i); r = rv(i);
    if r==0, Jc{i+1} = zeros(P.mres,0); continue; end
    Jc{i+1} = 2*(P.W*(P.Ssym(:,P.rows{i})*kron(sparse(Y{i}),speye(N))));
end
J = [Jc{:}];
end
