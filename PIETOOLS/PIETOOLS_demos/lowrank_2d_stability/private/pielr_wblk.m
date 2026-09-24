function v = pielr_wblk(w,Ns,rv,i,Kf)                                       % CC, 09/23/2026
% PIELR_WBLK  Slice block i of the factor out of the LM variable vector.
%
%   w = [ z (Kf free) ; Y_1(:) ; ... ; Y_B(:) ]
%
% THE Kf TERM IS THE POINT.  pielr_certify's local wblk starts its offset at
% zero and sums Ns(j)*rv(j), i.e. it assumes w begins with Y_1.  bm_resid,
% which defines the layout, starts at k = P.Kf and reads the free coordinates
% out of w(1:Kf) first.  The two agree only when Kf == 0 -- true for every
% stability program, and false for l2gain, where gamma and the free lpivar
% operator give Kf = 44 on the transport example.  With Kf > 0 the local wblk
% returns a window shifted by Kf entries into the wrong place, silently, and
% the face built from it is a face of the wrong subspace.

if nargin<5 || isempty(Kf), Kf = 0; end
k = Kf;
for j = 1:i-1, k = k + Ns(j)*rv(j); end
v = w(k+(1:Ns(i)*rv(i)));
end
