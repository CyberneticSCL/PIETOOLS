function w = pielr_fitw(Y0,rv,Ns,z0)                                        % CC, 09/23/2026
% PIELR_FITW  Assemble the LM variable vector w = [ z ; Y_1(:) ; ... ; Y_B(:) ]
% from a per-block factor, fitting it to the current rank profile.
%
% Generalises pielr_certify's local fitw in two ways.
%
% 1. IT CARRIES THE FREE BLOCK z.  fitw returns the Y blocks alone, which is
%    correct only when Kf == 0.  bm_resid reads w(1:Kf) as the free
%    coordinates, so with Kf > 0 (l2gain: gamma plus the free lpivar operator)
%    a w without them is misaligned by Kf entries from its first element on.
%
% 2. IT IS USED FOR THE RANK-LADDER WARM START, not only for opts.w0.  The
%    surplus columns are filled with SMALL RANDOM values, never zeros: a zero
%    column of Y_i has an identically zero Jacobian column block, because that
%    block is 2*W*Ssym(:,rows_i)*kron(Y_i(:,c),I) and kron(0,I) = 0.  The LM
%    step therefore has exactly no component along those coordinates, in both
%    the primal and the dual branch, and the added rank is inert for the whole
%    run -- the rung reproduces the previous rung's point.  pielr_certify's
%    padw pads with zeros and sits on precisely that fixed point; the
%    package's own 1-D test T4 pins it as defect B2 ("zero pad max||J(:,padded)||
%    = 0, EXPECTED EXACTLY 0") and shows this fill is mobile.  padw is left
%    alone there so pielr_certify remains the behaviour baseline.

if nargin<4, z0 = []; end
w = z0(:);
for i = 1:numel(Ns)
    Yi = zeros(Ns(i),rv(i));
    Y  = Y0{i};
    kk = min(size(Y,2),rv(i));
    if kk>0, Yi(:,1:kk) = Y(:,1:kk); end
    if rv(i) > kk
        s = norm(Y,'fro')/max(sqrt(numel(Y)),1);       % rms entry of the seed
        if ~isfinite(s) || s<=0, s = 1; end
        Yi(:,kk+1:end) = 1e-3*s*randn(Ns(i),rv(i)-kk);
    end
    w = [w;Yi(:)]; %#ok<AGROW>
end
end
