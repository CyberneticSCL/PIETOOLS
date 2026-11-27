function [q,V,dA] = bm_proj(P,q,rv,mode)
% mode 'affine': orthogonal (Frobenius) projection onto {X : A(X)=b}.
%   F = W*(A(q)-b) and  A^+(A(q)-b) = Ssym'*(W'*F), with ||A^+(...)|| = ||F||,
%   so the affine projection is exact and costs two sparse matvecs.
% mode 'rank'  : projection onto {X >= 0, rank(X_i) <= rv(i)} block by block.
switch mode
case 'affine'
    F = P.W*(P.Ssym*q - P.bf);
    q = q - P.Ssym'*(P.W'*F);
    dA = norm(F);  V = [];
case 'rank'
    B = numel(P.Ns); V = cell(1,B);
    for i=1:B
        N=P.Ns(i); Q=reshape(q(P.rows{i}),N,N); Q=(Q+Q')/2;
        [Wv,Dv]=eig(Q); dv=diag(Dv);
        [dv,ord]=sort(dv,'descend'); Wv=Wv(:,ord);
        k=min(rv(i),N); kp=false(N,1); kp(1:k)=dv(1:k)>0;
        V{i}=Wv(:,kp)*diag(sqrt(dv(kp)));
        Qn=V{i}*V{i}'; q(P.rows{i})=Qn(:);
    end
    dA = norm(P.W*(P.Ssym*q - P.bf));
end
end
