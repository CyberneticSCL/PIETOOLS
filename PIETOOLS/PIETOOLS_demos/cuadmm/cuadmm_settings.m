function S = cuadmm_settings(tier,dim)                                      % CC, 09/25/2026
% CUADMM_SETTINGS  A tiered PIETOOLS settings ladder for the cuADMM (GPU,
% first-order ADMM) arm.  Same signature and fields as pielr_settings, so a
% caller can swap one for the other.
%
%   S = cuadmm_settings(2)       % 1-D, tier 2
%   S = cuadmm_settings(2,2)     % 2-D
%
% STATUS: WRITTEN, NOT YET VALIDATED.  Every tier below follows pielr_settings,
% plus cuADMM-specific fields grounded in the 2026-09-24 baseline.  No tier has
% been run end to end through this file.
%
% Initial coding CC - 09/25/2026, after consulting the low-rank session on
% which of pielr_settings' findings are properties of the LPI (and so transfer)
% and which are properties of Burer-Monteiro (and so do not).
% CC, 09/26/2026: psatz paragraph rewritten from "necessity not isolated" to
%   the low-rank session's preset x flag cross: necessary for a searched P in
%   1-D, inert at 'stripped' (1x1 psatz block). Tier table unchanged -- it
%   already had psatz off at tier 0 and on from 'light' up.
% CC, 09/26/2026: S.cuadmm.bisect validated_to_m 273 -> 0. Its one validation
%   was h2c_rd1 (m = 273; H2, dropped from the regime today). The banked 1-D
%   Hinf rungs refute it at m = 243 and 268: at 0.99 gam* pinf is 3.8e-4
%   (hinf_rd1) and 3.2e-4 (hinf_rd1_hv), under 1e-3, where Mosek returns
%   PRIMAL_INFEASIBLE_CER. It agrees on hinfdu_rd1 (1.5e-3 vs 4.5e-4).
%
% ---------------------------------------------------------------------------
% WHAT TRANSFERS FROM pielr_settings, AND WHY
%
%  TIER BOUNDARIES -- transfer unchanged.  pielr_settings is an ACCURACY ladder,
%  not a cost ladder: even veryheavy is a 156x318 Jacobian in 1-D, trivially
%  affordable for BM, so its tier boundaries sit where the accuracy changes.
%  Accuracy is a property of the program, not the solver.
%
%  eppos = eppos2 = 1e-2, PINNED ON EVERY TIER -- transfers, and should matter
%  MORE here.  At eppos ~1e-6 the trivial point X=0 satisfies the equality rows
%  to ~1e-9 and acceptance sits on the noise floor; at 1e-2 that pathology
%  disappears, and the gate's ratio metric is invariant across six orders of
%  eppos, so raising it does not flatter the residual (low-rank session,
%  measured).  cuADMM's stopping rule divides by (1+||b||), which makes the
%  scale of b load-bearing: the only baseline case with an O(1) right-hand side
%  (wellposed_rd1, ||b|| = 17.7) converged in 44 iterations, against hundreds
%  to thousands for the eppos2 ~ 1e-6 feasibility cases.
%
%  THE PSATZ FLAG -- ON in tiers 1-3.  In 1-D it is NECESSARY, and it is INERT
%  below 'light'.  Measured by the low-rank session (2026-09-25/26) with the
%  full preset x flag cross, a searched P, Dup 1, eppos 1e-2 and solver='best'
%  (every solver on the path, the smaller residual kept); not re-run here:
%
%      lam/lam*   stripped OFF  stripped ON   light OFF   light ON
%       0.10       FAILS         FAILS         FAILS       certifies
%       0.50       FAILS         FAILS         FAILS       certifies
%       0.90       FAILS         FAILS         FAILS       FAILS
%
%    - light+OFF fails at every lam, so the psatz is necessary for a searched P,
%      not only for P = I (where poinc1d.m showed the same thing: reach 0.00
%      off at Dup 1-3 versus 0.50/0.9999 on).
%    - stripped+ON is BIT-IDENTICAL to stripped+OFF.  At 'stripped' the psatz
%      block is 1x1 ([6 8 1] against light's [10 15 6]): a multiplier with no
%      monomials cannot correct anything, so the flag is a no-op there.  Preset
%      and flag INTERACT: a tier that turns the psatz on at too coarse a base
%      pays for the extra block, gains nothing, and looks like evidence the
%      psatz does not help.  Do not enable it below 'light' -- hence tier 0 has
%      it off.
%    - Every FAILS cell above returns the trivial point X = 0: the residuals are
%      identical to 5 digits across presets and solvers.  So they read "no
%      certificate found", not "a near-miss".
%    - 2-D differs: psatz-off with a searched P certifies to ~0.25*lam* under
%      SeDuMi (7.3e-06 at 0.10*lam*), and the linear generators EXTEND reach
%      (1e-08 to 0.90*lam*) rather than enable it.  And in 2-D this ladder's
%      flag does not reach the executive at all -- see the 2-D note below.
%
%  Dup ON dd2/dd3 AS THE ACCURACY LEVER -- carried over as the ladder's axis.
%
% WHAT DOES NOT TRANSFER
%
%  lmit_min is a BM artefact (LM iterations).  Its analogue here is the cuADMM
%  iteration budget below, calibrated on cuADMM's own data, not BM's.
%
% ---------------------------------------------------------------------------
% WHAT IS DIFFERENT ABOUT cuADMM, measured 2026-09-24
%
%  1. ITERATION COUNT, NOT SIZE, IS THE COST.  t/iter ~ m^1.10 and flat in
%     sum N^3; dense eig <= 6% of time.  Memory is nearly free: < 1 GB of GPU
%     at m = 138,240, where Mosek's dense Schur complement alone needs 142 GB.
%  2. nnz(At) COSTS MORE THAN m.  Each iteration is dominated by the SpMV, so
%     the 2-D linear psatz generators, which add +30% to m but 6.6x to nnz(At)
%     (557k -> 3.66M), cost cuADMM roughly in proportion to the nnz.  Hence
%     they are an OPTION below, not a tier default.
%  3. OBJECTIVE FORM IS GAP-LIMITED.  1 of 18 objective cases converged at
%     1e-4, against 19 of 19 feasibility cases, typically already feasible to
%     ~1e-5 while the duality gap stalls.  So this ladder is defined on the
%     FEASIBILITY form; objective classes are run by bisection on a pinned
%     gamma (bl_fix.m appends e_j'x = gam and zeros c).
%  4. A HEAVIER BASE PRESET CAN CONVERGE IN FEWER ITERATIONS.  hinf_rd1_hv
%     (heavy) converged at 1e-4 in 4,184 iterations where hinf_rd1 (light) was
%     still gap-limited past 10,300.  light and heavy differ ONLY in n_order1,
%     i.e. on the BASE-PRESET axis, not on Dup.  The low-rank arm sees the same
%     split: raising the base preset at fixed Dup sometimes REDUCES time (tier 3
%     beat tier 2 on 7 of 18 cases), while raising Dup reliably increases it
%     (tier 1 fastest on 18 of 18).  NOT built into the tiers: two data points
%     on one class, and it may be a property of cuADMM's convergence rather
%     than of the program.  Check before skewing the ladder heavier.
%
% ---------------------------------------------------------------------------
% THE LADDER
%
%   tier  base      Dup  psatz  cuADMM tol  iter budget  intent
%    0    stripped   1    0       1e-4        20,000     smoke / reachability
%    1    light      1    1       1e-4        20,000     default
%    2    light      2    1       1e-6       100,000     accuracy
%    3    heavy      2    1       1e-6       100,000     accuracy, heavier base
%
%  ITERATION BUDGETS are set from the largest counts observed on FEASIBILITY
%  cases in the baseline, with headroom: at 1e-4 the maximum was 4,574
%  (nonlinear Fisher, m ~ 5000); at 1e-6 it was 33,236 (stabpde_rd1).  The
%  sample is the 19 feasibility cases up to m = 138,240; a larger or worse-
%  conditioned program may need more.  Objective cases did NOT inform these
%  numbers -- they need bisection, not a bigger budget.
%
%  cuADMM's own options (sig, sigscale, sig_thresh, switch_admm) are at STOCK
%  values.  They have never been swept on PIETOOLS programs.
%
% IN 2-D THIS LADDER COLLAPSES -- inferred from the code, not measured.  The
% 2-D executives REPLACE the whole settings struct with settings_2d and keep
% only sos_opts (PIETOOLS_stability_2D.m:64: settings = settings.settings_2d).
% So the tier's 1-D fields -- dd2, dd3, override2 -- are discarded in 2-D, and
% only the base preset's settings_2d and the eppos written into it survive:
%    tiers 1 and 2 pose the IDENTICAL 2-D program (both light; Dup never
%                  reaches settings_2d), and
%    the psatz flag is INERT in 2-D: eq_use_psatz stays at the shipped [0;0]
%                  in every tier. Use cuadmm_linear_psatz to set it.
% pielr_settings has the same property.  Giving the tiers meaning in 2-D needs
% 2-D degree bumps on settings_2d (LF_deg, eq_deg), which have not been
% measured and are not invented here.
%
% THE BISECTION DECISION RULE IS THE HIGHER-PRIORITY OPEN ITEM.  cuADMM has no
% infeasibility certificate and on pinned-gamma programs neither side converges,
% so feasibility is read from the primal residual PLATEAU at a fixed budget.
% That rule (pinf < 1e-3 after 5,000 iterations) reproduced Mosek's certificates
% on h2c_rd1 and hinfdu_rd1, and FAILED on hinf_rd1 and hinf_rd1_hv (0.99 gam*
% called feasible, m = 243 and 268) and on a control at m = 14,006: it called
% gamma = 0.70x the analytic gain feasible, because pinf decays smoothly through
% the true boundary there.  A short budget gives a loose answer; a broken rule
% gives a WRONG one.  S.cuadmm.bisect records that it is validated at no m; it
% must not drive a bisection without a stagnation test and a known-answer
% control.  (CC, 09/26/2026: paragraph corrected, see header log.)
%
% INPUT
%   tier   0,1,2,3 (default 1)
%   dim    1 or 2 (default 1)
% OUTPUT
%   S      lpisettings struct with the tier applied, plus S.cuadmm (solver
%          options, budget, form, bisection rule), S.poincare, and the
%          provenance fields S.tier, S.tier_base, S.tier_Dup, S.dim.

if nargin<1 || isempty(tier), tier = 1; end
if nargin<2 || isempty(dim),  dim  = 1; end
assert(ismember(tier,0:3),'cuadmm_settings: tier must be 0, 1, 2 or 3');
assert(ismember(dim,1:2),'cuadmm_settings: dim must be 1 or 2');

% ---- tier table: identical to pielr_settings (accuracy boundaries transfer)
switch tier
    case 0, base = 'stripped';  Dup = 1;  psatz = 0;
    case 1, base = 'light';     Dup = 1;  psatz = 1;
    case 2, base = 'light';     Dup = 2;  psatz = 1;
    case 3, base = 'heavy';     Dup = 2;  psatz = 1;
end
S = lpisettings(base);

% sossolve would otherwise overwrite At/b/K/c with the sospsimplify-reduced
% system (sossolve.m:256-262), so the dumped SDP would not be the one solved.
S.sos_opts.simplify = false;
% The reference solve is SeDuMi, as in pielr_settings.  Mosek is faster but was
% measured to FAIL on 2-D psatz-off stability (rel_b 1.5e-03 / 2.18 / 2.08 at
% 0.10 / 0.25 / 0.50 lam*) where SeDuMi certifies (7.3e-06 at 0.10 lam*).
S.sos_opts.solver = 'sedumi';

% ---- eppos pinned at 1e-2 on every tier (see header)
S.eppos = 1e-2;   S.eppos2 = 1e-2;
if isfield(S,'settings_2d'), S.settings_2d.eppos = 1e-2*ones(4,1); end

% ---- negativity psatz and the Dup bump on the slack -- as pielr_settings
S.override2 = double(~psatz);          % 1 means OFF: the flag is inverted
n1 = 1;  n2 = 1;  n3 = 1;              % the 1-D light/heavy base orders
if strcmp(base,'heavy'), n1 = 2; end   % light and heavy differ ONLY here
n4 = n2 + n3;
S.dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup],   [n2+Dup-1, n3+Dup,   n4+Dup]};
S.dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};
S.poincare = struct('psatz',psatz,'dbump',max(0,tier-1));

% ---- 2-D linear psatz generators: an OPTION, off by default in every tier;
% enable with S = cuadmm_linear_psatz(S). Not stored inside settings_2d, because
% the 2-D executives REPLACE the whole settings struct with settings_2d
% (PIETOOLS_stability_2D.m:64) and an extra field would ride along into it.

% ---- cuADMM solver options
tolv   = [1e-4  1e-4  1e-6   1e-6];
budget = [20000 20000 100000 100000];
S.cuadmm = struct( ...
    'tol',         tolv(tier+1), ...
    'max_iter',    budget(tier+1), ...
    'eig_streams', 15, ...          % stock
    'sig',         1e2, ...         % stock, never swept
    'sigscale',    2.0, ...         % stock, never swept
    'sig_thresh',  500, ...         % stock, never swept
    'switch_admm', 0, ...           % stock: never switch to sGS-ADMM
    'form',        'feasibility', ...
    'linear_psatz', false, ...      % 2-D only; see cuadmm_linear_psatz
    'ld_library_path', '/usr/lib/wsl/lib');   % MANDATORY on the WSL build
% The bisection rule, with the range it has been validated over: none. It
% misclassifies at m = 243, 268 and 14,006 (see header); do not use it.
% CC, 09/26/2026: was validated_to_m 273, refuted_at_m 14006, note 'needs a
% stagnation test and a known-answer control above validated_to_m'.
S.cuadmm.bisect = struct( ...
    'rule',           'pinf_plateau', ...
    'budget',         5000, ...
    'pinf_threshold', 1e-3, ...
    'validated_to_m', 0, ...
    'refuted_at_m',   [243 268 14006], ...
    'note', 'validated at no m; needs a stagnation test and a known-answer control'); % CC, 09/26/2026

S.tier = tier;   S.tier_base = base;   S.tier_Dup = Dup;   S.dim = dim;
S.ladder = 'cuadmm';
end


