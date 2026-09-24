function R = pielr_opcheck(prog,H,P,q,A)                                    % CC, 09/23/2026
% PIELR_OPCHECK  Operator-level verification of an LPI candidate, in either
% spatial dimension and for any adapted executive.  THE acceptance test.
%
% Pushes the Gram vector q back through the real PI operators, rebuilds the
% executive's own constraint expression from the substituted operators, and
% measures the residual against a normaliser the adapter names.
%
% WHY THE NORMALISER CHANGED (the one behaviour change in this rewrite).
% opcheck_2d reports rel = max|Res| / max|Qop|, dividing by the operator that
% must VANISH.  tests_1d/opcheck.m divides by max|Pop| instead and states why:
%     "rel is INVALID wherever the candidate drives Dop -> 0 as an operator:
%      it then reads 0/0 and improves for a purely artifactual reason
%      (measured on the Example 21 beam, max|Dop| = 2.4e-10 manufactured a
%      clean degree threshold that does not exist).  relP divides by max|Pop|
%      instead, which cannot collapse because Pop >= eppos2*I by
%      construction.  Accept on relP, never on rel alone."
% So the package's two halves accepted on different quantities, and the 2-D
% half accepted on the one its own 1-D half documents as invalid.  A benchmark
% spanning both dimensions cannot be built on that.  Here the adapter supplies
% the normaliser together with the reason it is bounded below (A.resid's
% S.Nrm / S.why), R.rel is measured against it, and the old 2-D quantity is
% still reported as R.rel_d so existing numbers stay comparable.
%
% WHAT IS *NOT* CHANGED.  The second half of the gate -- every Gram block PSD
% -- and its tolerance are exactly opcheck_2d's.  Worth knowing when reading a
% pass: on every acceptance path this package ships, the candidate is either
% Y*Y' or a V*S*V' with S projected onto the PSD cone, so this half is
% satisfied by construction and cannot fail; the whole verdict rests on R.rel.
%
% THE THRESHOLD IS RELATIVE TO A REFERENCE (CC, 09/23/2026).  The acceptance
% threshold used to be the constant 1e-6.  MEASURED, on the historical 2-D
% benchmark (nb_rd2d(1,0.50*lam*), set2d_deg(4,[]), N=[8 744], m=5920): the
% INTERIOR-POINT solution of that program -- the full-rank answer -- reaches
% rel = 3.497e-06.  A 1e-6 threshold there demands more precision than the
% program admits from any method, so it is not an acceptance criterion but a
% guarantee of failure.  Across tier 1 the interior-point residual spans
% 4.15e-10 to 9.18e-07, three and a half orders of magnitude, which no single
% absolute constant can span.
%
% So the gate is now
%       R.score = R.rel / ref                    <- the number to report
%       R.ok    = R.rel <= R.thresh  AND  every Gram block PSD
%       R.thresh = max(abs, k*ref)   when ref is credible (ref <= refmax)
%                = abs               otherwise
% where ref is the residual a reference solve achieves on the SAME program.
% The 'abs' floor is kept at 1e-6 by default so that nothing which passed
% before fails now: every certificate this package has accepted has
% rel < 1e-6 <= max(1e-6, k*ref).  Setting abs = 0 gives the pure ratio gate.
%
% WHY A FLOOR AT ALL, rather than the ratio alone.  The reference residual
% measures how well the reference solver converged, not a floor on achievable
% accuracy: on easy programs it reaches 4e-10 and a pure ratio gate would then
% demand 4e-10*k of the low-rank search, i.e. it would be HARSHEST exactly
% where the problem is easiest.  MEASURED ratios of the currently-accepted
% certificates run from 6.2 to 3908, so a pure ratio gate either accepts
% everything (k ~ 1e4) or rejects almost everything (k ~ 1e2).  The floor is
% what makes the relative part a relaxation on hard programs instead of a
% tightening on easy ones.
%
% CHOOSING k, MEASURED (tier 1, closed loop -- k changes which point the search
% RETURNS, since the ladder stops at the first rung that certifies and the
% gamma bisection keeps whatever certifies, so this cannot be answered by
% re-scoring banked points):
%     k      certs  binds  score med   gamma ratios
%     1      11     0      81.29       1.0819 1.1091 1.0010 2.8264
%     10     11     2      81.29       1.0819 1.1091 1.0010 1.6856
%     100    13     6      95.39       1.0779 1.1083 1.0010 0.6143
%     1000   14     11     667.6       1.0148 1.0410 1.0010 0.1936
% The gamma column FALSIFIES the large values.  A face restriction can only
% LOSE feasible points, so the restricted optimum cannot be below the full
% optimum and a gamma ratio under 1 is a contradiction, not tolerance.  At
% k = 100 and 1000 the gate is loose enough to accept points that do not
% satisfy the LPI and the bisection then drives gamma to 0.61x and 0.19x the
% reference.  k <= 10 is validated, k >= 100 is refuted, by the suite's own
% data.  The l2gain cases therefore act as a standing calibration of k, which
% is a reason to keep at least one of them in any run that changes the gate.
%
% With no reference supplied the behaviour is exactly the old absolute gate.
%
% UNITS.  q is in NORMALISED-b units, as produced by bm_setup / bm_resid; this
% routine multiplies by P.nb0 internally.  Passing original units silently
% reports rel = 1.
%
% INPUT   prog  the assembled program (solved or not; only solinfo is used)
%         H     operator handles from A.build
%         P     bm_setup package (needs .nb0, .rows; .N or .Ns)
%         q     length P.Ntot vector, normalised-b units
%         A     adapter from pielr_lpi
% OUTPUT  R.rel      max|Res| / max|Nrm|          <-- THE number
%         R.rel_d    max|Res| / max|Res-operator|  (opcheck_2d's old quantity)
%         R.maxRes .maxNrm .maxDop  the maxima the ratios are built from
%         R.nrm_why  why the denominator cannot collapse
%         R.mineig(i) .normQ(i) .rank(i)  per Gram block, ORIGINAL units
%         R.psd      true iff every block's min eig >= -1e-8 * its own norm
%         R.score    R.rel / ref, or NaN when no reference was supplied
%         R.thresh   the threshold actually applied
%         R.ok       true iff R.rel <= R.thresh AND R.psd
%         R.aux      adapter extras (gamma, the rebuilt operators)
%
% A.gate, when the adapter carries one: .abs (absolute floor, default 1e-6),
% .k (multiple of the reference, default 10), .ref (the reference residual,
% NaN when none).  pielr_solve fills it; pielr_lpi does not, because the
% reference is a property of the PROGRAM, not of the executive.

pg = prog;
pg.solinfo.info = struct('verified_by','pielr_opcheck');
pg.solinfo.RRx  = full(q(:))*P.nb0;        % back to the original b scaling

S = A.resid(pg,H);

% THE GATE MUST COVER EVERY OPERATOR EQUALITY THE PROGRAM IMPOSES, not only
% the one the analysis is about (CC, 09/23/2026).  An adapter may therefore
% return S.Res / S.Nrm as CELLS, one pair per equality, and the verdict is the
% WORST of the per-equality ratios.
%
% WHY.  Raised by the parallel cuADMM session and checked here against every
% builder in this package: build_l2gain_1d imposes TWO equalities -- the
% coupling Top'*Qop = Rop and the negativity Deop + Dop = 0 -- and l2g_resid
% rebuilt only the second.  A point can then satisfy the negativity identity
% while Rop is not Top'*Qop at all, i.e. the operator asserted positive
% semidefinite is not the operator appearing in the Lyapunov inequality, and
% the gate would accept it.  build_stab_1d, build_stab_2d_st2 and
% build_poincare_1d impose exactly one equality each and are unaffected.
%
% Their own case is sharper and does not apply here, but is worth recording:
% PIE2PDEstability puts the eppos constant in Pop, which appears ONLY in the
% coupling equality, so its negativity residual is identically zero at X = 0
% and NO choice of denominator on that equality can reject the trivial point.
% This package's 'stability' adapter builds the DIRECT form, where
% Dop(0) = T'(eppos I)A + A'(eppos I)T is nonzero, so it is not exposed --
% but anything routed through the Q form would be.
Res = S.Res;   Nrm = S.Nrm;
if ~iscell(Res), Res = {Res};  end
if ~iscell(Nrm), Nrm = {Nrm};  end
assert(numel(Res)==numel(Nrm), ...
    'pielr_opcheck: %d residual(s) against %d normaliser(s)',numel(Res),numel(Nrm));
np_ = numel(Res);
R.maxRes_i = zeros(1,np_);  R.maxNrm_i = zeros(1,np_);  R.rel_i = zeros(1,np_);
for i = 1:np_
    R.maxRes_i(i) = mx_(Res{i});
    R.maxNrm_i(i) = mx_(Nrm{i});
    R.rel_i(i)    = R.maxRes_i(i)/max(R.maxNrm_i(i),realmin);
end
[R.rel,iw] = max(R.rel_i);          % the worst equality decides
R.maxRes  = R.maxRes_i(iw);
R.maxNrm  = R.maxNrm_i(iw);
R.worst   = iw;
R.nrm_why = S.why;
% the quantity opcheck_2d accepted on, kept so old logs remain comparable and
% so a collapsing denominator is visible rather than silent
if isfield(S,'aux') && isfield(S.aux,'Dop')
    R.maxDop = pielr_maxop(S.aux.Dop);
else
    R.maxDop = R.maxNrm;
end
R.rel_d = R.maxRes/max(R.maxDop,realmin);
R.aux   = S.aux;

% ---- per-block spectrum of the Gram blocks (unchanged from opcheck_2d) ---
B = numel(P.rows);
R.normQ = zeros(1,B);  R.mineig = zeros(1,B);  R.rank = zeros(1,B);
for i = 1:B
    N = round(sqrt(numel(P.rows{i})));
    Q = reshape(q(P.rows{i}),N,N);
    Q = (Q+Q')/2;
    e = sort(eig(Q),'descend');
    R.normQ(i)  = norm(Q,'fro')*P.nb0;
    R.mineig(i) = e(end)*P.nb0;
    R.rank(i)   = sum(e > max(e(1),eps)*1e-9);
end
R.psd = all(R.mineig >= -1e-8*max(R.normQ,realmin));

% ---- the threshold ------------------------------------------------------
% THE REFERENCE MUST ITSELF BE CREDIBLE (CC, 09/23/2026).  A first cut applied
% max(abs, k*ref) unconditionally and accepted THREE garbage points, caught by
% the suite's negative controls: on programs the reference solver cannot solve
% either it returns ref = 8.883, 1.500, 3.521, so k*ref is 88.8, 15.0, 35.2
% and residuals of 4.043e-01, 6.285e-02, 7.180e-02 "passed".  Ten times
% garbage is garbage.  The relative branch exists to RELAX the threshold on a
% program that cannot reach 1e-6, not to abolish it on a program nobody can
% solve, so it only applies while the reference is itself in a credible range.
%
% g.refmax is where that line sits, and THE BINDING QUANTITY IS k*refmax, not
% k alone: that product is the loosest threshold this gate can ever apply, so
% it is what has to stay credible.
%
% refmax was first set to 1e-4 on tier-1 evidence -- credible references ran to
% 9.181e-07, non-credible ones started at 1.500, and 1e-4 sat in the gap.  THE
% 2-D CASES POPULATED THAT GAP and the choice failed: rd2d-deg3 at frac 0.50 is
% a NEGATIVE CONTROL whose reference reaches only 6.283e-05, and with
% refmax = 1e-4 the relative branch applied, the threshold became 6.283e-04,
% and a point at rel = 1.020e-04 was accepted.  An operator identity holding to
% four relative digits is not a certificate, whatever its ratio to the
% reference (its score was 1.623, the BEST in the suite -- which is the lesson:
% a good ratio on a program nobody can solve is not a good result).
%
% refmax = 1e-5 is measured from both sides:
%   admit  2-D deg 4 frac 0.50, ref 3.497e-06   <- the case that motivated all
%   admit  2-D deg 3 frac 0.10, ref 1.900e-06      of this
%   REJECT 2-D deg 3 frac 0.50, ref 6.283e-05   <- the negative control
% and it caps the loosest possible threshold at k*refmax = 1e-4 for k = 10.
g = struct('abs',1e-6,'k',10,'ref',NaN,'refmax',1e-5);
if isfield(A,'gate') && ~isempty(A.gate)
    fn = fieldnames(A.gate);
    for i = 1:numel(fn), g.(fn{i}) = A.gate.(fn{i}); end
end
if isfinite(g.ref) && g.ref > 0
    R.score = R.rel/g.ref;       % the score is reported whatever the threshold
else
    R.score = NaN;
end
if isfinite(g.ref) && g.ref > 0 && g.ref <= g.refmax
    R.thresh = max(g.abs, g.k*g.ref);
else
    R.thresh = g.abs;            % no usable reference: the absolute gate
end
R.gate_ref = g.ref;   R.gate_k = g.k;   R.gate_abs = g.abs;
R.gate_refmax = g.refmax;
R.ok = (R.rel <= R.thresh) && R.psd;
end

% =========================================================================
function m = mx_(X)                                                         % CC, 09/23/2026
% max|coefficient| of an operator, or a plain nonnegative scalar passed
% straight through.  A normaliser is sometimes a number rather than an
% operator -- the l2gain coupling uses max(|Top'Qop|,|Rop|) -- and
% manufacturing an operator just to carry it would be worse.
if isnumeric(X) && isscalar(X), m = abs(double(X)); return, end
m = pielr_maxop(X);
end
