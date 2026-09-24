function R = bm_report(w,P,rv)
% Diagnostics for a Burer-Monteiro point, in BOTH metrics.
[F,~,q] = bm_resid(w,P,rv);
raw = P.Ssym*q - P.bf;
R.dist_affine = norm(F);                  % Frobenius distance of Q to {A(X)=b}
R.raw_abs     = norm(raw);                % ||A(YY') - b||
R.raw_rel     = norm(raw)/P.nb;           % ||A(YY') - b|| / ||b||     <- the requested metric
R.raw_relinf  = max(abs(raw))/max(abs(P.bf));
B = numel(P.Ns);  R.normQ = zeros(1,B);  R.mineig = zeros(1,B);  R.numrank = zeros(1,B);
for i=1:B
    Q = reshape(q(P.rows{i}),P.Ns(i),P.Ns(i)); Q=(Q+Q')/2;
    e = sort(eig(Q),'descend');
    R.normQ(i)=norm(Q,'fro'); R.mineig(i)=e(end);
    R.numrank(i) = sum(e > max(e(1),eps)*1e-9);
end
R.normQtot = norm(R.normQ);
end
