function E = pielr_eta(At,b,x,P)                                            % CC, 09/27/2026
% PIELR_ETA  Row-normwise backward error of a decision vector, and of its
% PSD-clipped form.  REPORTED ALONGSIDE the existing residuals; it decides
% nothing yet.
%
% WHY NOT ||At'x-b||/||b||.  Measured on this package's own programs: b is
% 86.4-98.8% ZEROS -- 2 nonzeros in 130 rows on gain-transport, 2 in 160 on
% gain-reacdiff, 70 in 3456 on rd2d-deg3-f010.  A ratio of 2-norms is then an
% ABSOLUTE test dominated by the handful of rows carrying b, it is not
% invariant to row scaling, and it tightens like 1/sqrt(m) as the program
% grows.  Row scaling is not uniform here either: ||A_i||_2 spans 480x on
% rd1d-lam0.5, 561x on gain-reacdiff and 1120x on advdiff1d, so the same point
% passes or fails depending on how the rows happen to be scaled.
%
% THE REPLACEMENT.
%       eta = max_i |r_i| / (||A_i||_2 * ||x||_inf + |b_i|),   r = At'x - b
% Invariant to row scaling, and it means "x solves EXACTLY a problem in which
% every row's data is perturbed by relative <= eta".  Pairing ||A_i||_2 with
% ||x||_inf overstates the true backward error, so the number is conservative.
%
% WHY THE GLOBAL ||x||_inf AND NOT A PER-ROW SCALE.  The componentwise
% Oettli-Prager omega_i = |r_i|/((|A_i| |x|)_i + |b_i|) is exactly 1 on the
% forced-zero Gram rows (b_i = 0 with same-sign diagonals): there the exact
% solution is 0 and any solver returns ~1e-8, so omega reports 1 for a point
% that is fine.  That is the missing strict interior showing up, not a signal.
%
% THE PSD CLIP.  A Gram block is only a certificate if it is in the cone, so
% clip each block to its PSD part, X+ = V*max(L,0)*V', and report eta(X+).
% The tolerated negative eigenvalue converts roughly 1:1 into row error, so a
% point reported at eta ~1e-17 with lambda_min/lambda_max ~ -3e-9 has
% eta(X+) ~4e-9 -- nine orders larger than the headline number.  eta(X+) is
% the backward error of the WHOLE certificate and is the honest one.
%
% ETA ALONE IS UNSOUND -- do not gate on it without a scale guard.  It is
% relative to ||x||_inf, so a point growing along a near-recession direction
% flatters it: the cuADMM session measured SeDuMi infeasible points at
% 0.9966-0.9982 gamma_F with ||X|| ~2e3 against ~10 for genuine points, and
% eta(X+) 8.5e-08 to 1.2e-07 -- they would PASS a 1e-7 tolerance.  E.normx is
% returned so a caller can apply a norm guard relative to a certified
% reference point; an absolute eigenvalue floor also works but rejected 25
% genuinely feasible large-norm iterates in their set.
%
% WHAT IT DOES NOT FIX.  None of these programs has a strictly feasible point
% (a margin program max t s.t. Y_k >= t*s_k*I returned t* ~ 0 in every class,
% and in 2-D for every N=424 block).  So ANY residual test certifies a NEARBY
% problem, and separating feasible from infeasible near the boundary remains a
% calibration rather than a guarantee.  That is a property of the formulation,
% not of the metric.
%
% INPUT   At  Ntot x m  (constraint is At'*x = b),  b  m x 1
%         x   Ntot x 1  decision vector, ORIGINAL b scaling (not q)
%         P   the bm_setup package, for P.rows -- the per-block index sets
% OUTPUT  E.eta      row-normwise backward error of x
%         E.eta_psd  the same for the PSD-clipped x
%         E.normx    ||x||_inf, for the scale guard a gate would need
%         E.mineig   min eigenvalue over blocks, original units
%         E.clipnorm ||x - x_clipped||_inf, how much the clip moved the point
%
% COST: one mat-vec plus precomputed row norms, O(nnz(At)); the clip reuses an
% eigendecomposition the caller may already have done.

x = full(x(:));
E = struct('eta',NaN,'eta_psd',NaN,'normx',norm(x,inf), ...
           'mineig',NaN,'clipnorm',NaN);
rn = full(sqrt(sum(At.^2,1)))';          % ||A_i||_2, one per equality row
bf = full(b(:));
E.eta = etaof(At,bf,x,rn);

% ---- clip every Gram block to its PSD part ------------------------------
% The blocks are FULL N^2 column-major vectorisations indexed by P.rows{i},
% NOT packed triangles -- mirrored from pielr_opcheck rather than assumed, and
% an assumed svec packing here silently produced a meaningless eta_psd before
% this was checked.
xc = x;  me = inf;
for i = 1:numel(P.rows)
    idx = P.rows{i};
    n   = round(sqrt(numel(idx)));
    M   = reshape(xc(idx),n,n);   M = (M+M')/2;
    [V,L] = eig(M);   d = diag(L);
    me = min(me,min(d));
    if any(d<0)
        xc(idx) = reshape(V*diag(max(d,0))*V',[],1);
    end
end
E.mineig   = me;
E.clipnorm = norm(x-xc,inf);
E.eta_psd  = etaof(At,bf,xc,rn);
end

function e = etaof(At,b,x,rn)
r = At.'*x - b;
den = rn*norm(x,inf) + abs(b);
g = den > 0;                             % a zero row with b_i = 0 is vacuous
if ~any(g), e = 0; return, end
e = max(abs(full(r(g)))./full(den(g)));
end

