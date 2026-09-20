function P = bm_setup(Atf,bf,Ns,Kf,pre)
% Package the SeDuMi data for Burer-Monteiro.
%   x = [ z (Kf free) ; vec(Q_1) ; ... ; vec(Q_B) ],  Q_i = Y_i*Y_i'
%   constraint  Atf'*x = bf   (m equations)
%
% Ssym is Atf' with each block's (i,j)/(j,i) columns averaged, so that
% Ssym*vec(Q) = Atf'*vec(Q) for symmetric Q and the Jacobian simplifies to
%   J_i = 2 * Ssym(:,rows_i) * kron(Y_i, I_N).
%
% PRECONDITIONER (pre=1, default).  The raw residual is hopelessly scaled here:
% ||b|| ~ 5e-6 while feasible Q have ||Q||_F ~ 10, so A(Q) must cancel to 6
% digits before the residual is even O(||b||).  Rows of Ssym are symmetric
% matrices, so G = Ssym*Ssym' is the Gram of the constraint functionals in the
% Frobenius metric.  With G = U D U' and W = D^{-1/2}U' (null directions
% dropped), ||W*(Ssym*vec(Q)-b)|| is EXACTLY the Frobenius distance from Q to
% the affine set {X : A(X)=b}.  Same zero set, natural units.
% b is rescaled to UNIT NORM.  The LPI is nearly homogeneous (||b||~5e-6 comes
% only from the eppos*I term) so solutions span ||Q|| from ~1e-6 to ~10; with
% b normalised every residual reported is ALREADY ||A(Q)-b||/||b|| and absolute
% optimizer tolerances mean the same thing at every solution scale.
if nargin<5||isempty(pre), pre = 1; end
bf  = full(bf(:));
P.nb0 = norm(bf);          % original ||b||, for converting Q back
bf = bf/P.nb0;
S   = Atf';                       % m x Ntot
B   = numel(Ns);
rows = cell(1,B);  off = Kf;
for i=1:B
    rows{i} = off+(1:Ns(i)^2);  off = off+Ns(i)^2;
end
Ssym = S;
for i=1:B
    N = Ns(i);
    p = reshape(reshape(1:N^2,N,N)',[],1);        % vec(M') = vec(M)(p)
    Ssym(:,rows{i}) = 0.5*(S(:,rows{i}) + S(:,rows{i}(p)));
end
P.S = S;  P.Ssym = Ssym;  P.bf = bf;  P.Ns = Ns;  P.Kf = Kf;
P.rows = rows;  P.m = numel(bf);  P.Ntot = size(S,2);
P.nb = norm(bf);
P.pre = pre;
if pre
    G = full(Ssym*Ssym');
    G = (G+G')/2;
    [U,D] = eig(G);  d = diag(D);
    keep = d > max(d)*1e-12;
    P.W  = diag(1./sqrt(d(keep)))*U(:,keep)';     % k x m
    P.rankA = sum(keep);
    P.mres  = P.rankA;
    % b must lie in range(A); the dropped directions must carry no b-component
    P.bout = norm(bf - U(:,keep)*(U(:,keep)'*bf))/max(P.nb,eps);
else
    P.W = speye(P.m);  P.rankA = NaN;  P.mres = P.m;  P.bout = 0;
end
P.Wb = P.W*bf;
end
