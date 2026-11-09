function [V,best,hist] = bm_dr(P,rv,q0,maxit)
% Douglas-Rachford between  A = {X : A(X)=b}  and  B = {X >= 0, rank <= rv}.
%   y_{k+1} = y_k + P_A( 2 P_B(y_k) - y_k ) - P_B(y_k)
% Returns the best rank-rv PSD iterate seen and its distance to the affine set.
if nargin<4||isempty(maxit), maxit=400; end
y = q0; best = inf; V = []; hist = zeros(1,maxit);
for it=1:maxit
    [pb,Vb,dB] = bm_proj(P,y,rv,'rank');
    if dB < best, best = dB; V = Vb; end
    hist(it) = dB;
    pa = bm_proj(P,2*pb-y,rv,'affine');
    y  = y + pa - pb;
end
hist = hist(1:it);
end
