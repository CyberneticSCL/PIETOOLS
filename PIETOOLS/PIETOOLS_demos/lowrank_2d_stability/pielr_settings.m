function S = pielr_settings(tier,dim)                                      % CC, 09/23/2026
% PIELR_SETTINGS  A tiered PIETOOLS settings ladder sized for the low-rank
% (Burer-Monteiro) search.
%
%   S = pielr_settings(2)        % 1-D, tier 2
%   S = pielr_settings(2,2)      % 2-D
%
% WHY A NEW LADDER RATHER THAN THE STOCK PRESETS.  The stock presets vary the
% base degree n_order and little else, and they vary it on an axis that is
% nearly free for BM while leaving fixed the two things BM actually pays for.
% MEASURED, 1-D stability (x_t = x_ss + 0.5pi^2 x, Dirichlet, eppos 1e-2):
%
%   preset      ov2  blocks  N               m     Ntot   nv(r=2)
%   extreme      1     2     [2 4]            18      20     12
%   stripped     1     2     [6 8]            33     100     28
%   light        0     3     [10 15 6]        54     361     62
%   heavy        0     3     [11 20 11]       56     642     84
%   veryheavy    0     3     [36 71 52]      156    9041    318
%
% In 1-D every one of those is trivially affordable for BM: the Jacobian is
% m x sum(N_i r_i), so even veryheavy is 156 x 318.  The stock ladder is
% therefore not a cost ladder for BM in 1-D -- it is an ACCURACY ladder, and
% the accuracy lever is not where the preset names suggest.
%
% WHAT ACTUALLY MOVES ACCURACY, measured:
%
%  (1) PSATZ ON THE NEGATIVITY SLACK IS NECESSARY, not cosmetic.  With
%      override2 = 1 (psatz OFF -- note the inversion, 1 means off) the
%      program does not certify even at a tenth of the stability boundary:
%        stripped  lam/lam* 0.10  ipm rel 2.658e+00  FAILS
%        stripped  lam/lam* 0.50  ipm rel 5.290e+00  FAILS
%        light     lam/lam* 0.10  ipm rel 5.550e-09  certifies
%        light     lam/lam* 0.50  ipm rel 9.267e-08  certifies
%      Raising the plain degree cannot substitute for the boundary factor.
%      It also sets the BLOCK COUNT, which is what the rank ladder sweeps
%      across: 2 blocks without, 3 with.
%
%  (2) Dup ON dd2/dd3 IS THE ACCURACY LEVER, and it is cheap:
%        lam/lam* 0.50  Dup1 rel 9.267e-08 certifies | Dup2 4.111e-10 certifies
%        lam/lam* 0.90  Dup1 rel 8.883e+00 FAILS     | Dup2 4.599e-10 certifies
%        lam/lam* 0.99  Dup1 rel 9.771e+00 FAILS     | Dup2 4.866e-09 certifies
%      with m growing only 54 -> 59 (1.09x) and N [10 15 6] -> [10 26 15].
%      Reach goes from about half the boundary to essentially all of it for
%      9% more equality rows.
%
% (1) and (2) were reported to me by the parallel cuADMM session and are
% reproduced above on my own runs before being built on.
%
% THE LADDER.  Tier is an accuracy level; each tier is BM-affordable in 1-D,
% and the 2-D column is where the BM cost actually binds.
%
%   tier  negativity psatz  Dup  base      intent
%    0    OFF               1    stripped  deliberately too weak -- the
%                                          control that must give a LOOSE
%                                          bound or none at all
%    1    ON                1    light     the stock default
%    2    ON                2    light     the accuracy lever, ~9% more m
%    3    ON                2    heavy     both levers
%
% POINCARE is sized differently and gets its own sub-struct: that LPI comes
% through lpi_ineq, whose degrees are set by degbalance(P) from the operator
% rather than by any preset, so the only knobs are the psatz flag and a bump
% on the balanced degree.  S.poincare.psatz / S.poincare.dbump carry them.
%
% OUTPUT  S  a settings struct usable as pielr_solve's opts.settings, with
%            S.poincare attached for the Poincare adapter and S.tier recorded.

if nargin<1 || isempty(tier), tier = 1; end
if nargin<2 || isempty(dim),  dim  = 1; end
assert(ismember(tier,0:3),'pielr_settings: tier must be 0, 1, 2 or 3');

switch tier
    case 0, base = 'stripped';  Dup = 1;  psatz = 0;
    case 1, base = 'light';     Dup = 1;  psatz = 1;
    case 2, base = 'light';     Dup = 2;  psatz = 1;
    case 3, base = 'heavy';     Dup = 2;  psatz = 1;
end

S = lpisettings(base);
S.sos_opts.simplify = false;
% PIN THE SOLVER.  sossolve takes the first solver on the path and mosekopt
% leads its list, so an unset sos_opts.solver silently selects Mosek --
% verified present on this machine at C:\Program Files\Mosek\11.0.  A
% benchmark that does not pin this is not measuring what it says it is.
S.sos_opts.solver = 'sedumi';

% eppos: b is proportional to it for the stability LPI and nothing else
% contributes, so the stock 1e-4/1e-6 leaves the whole Lyapunov operator at
% epsilon scale.  Relative metrics are invariant to it (measured across six
% orders) but the Gram-level positivity margin is not.
S.eppos = 1e-2;   S.eppos2 = 1e-2;
if isfield(S,'settings_2d'), S.settings_2d.eppos = 1e-2*ones(4,1); end

% negativity psatz, and the degree bump on the slack
S.override2 = double(~psatz);          % 1 means OFF; see the header
n1 = 1;  n2 = 1;  n3 = 1;              % the 1-D light/heavy base orders
if strcmp(base,'heavy'), n1 = 2; end
n4 = n2 + n3;
S.dd2 = {n1+Dup,   [n2+Dup-1, n3+Dup,   n4+Dup],   [n2+Dup-1, n3+Dup,   n4+Dup]};
S.dd3 = {n1+Dup-1, [n2+Dup-2, n3+Dup-1, n4+Dup-1], [n2+Dup-2, n3+Dup-1, n4+Dup-1]};

% Poincare: psatz flag plus a bump on degbalance's own degrees
S.poincare = struct('psatz',psatz,'dbump',max(0,tier-1));

% BM BUDGET SCALES WITH THE TIER, and this is the number a caller most needs.
% VERIFIED on Poincare, whose answer is analytic (1/pi = 0.3183099), by
% sweeping lmit at each tier:
%
%   tier   lmit 400    lmit 2000   lmit 8000   interior-point on that tier
%    1     0.426909    0.426785    0.426785    0.4270904
%    2     0.342519    0.333951    0.318352    0.3183116
%    3     0.375286    0.318586    0.318392    0.3183102
%
% Three things to read off it.  The LADDER IS SOUND: the interior-point column
% tightens monotonically to 1/pi and tier 2 already reaches it to six digits.
% BM TRACKS THE LADDER once its budget is adequate -- 1/pi to four digits at
% tiers 2 and 3.  And the budget REQUIREMENT GROWS WITH THE TIER: tier 1 has
% converged by 2000 and does not move after, while tiers 2 and 3 are still
% cut at 400 (LM exits maxit:78 and maxit:73 there, falling to 16 and 14 by
% 8000).  Tier 3 at lmit 400 is WORSE than tier 2 at lmit 400 -- 0.375 against
% 0.343 -- purely because the larger blocks need proportionally more LM
% iterations.  Reading that as "tier 3 is worse" would be a truncation
% artefact, which is the same trap that produced this package's earlier reach
% figures.
%
% Note also the tier-1 ceiling is the SETTINGS' and not the algorithm's: BM
% and the interior-point solve agree there to 0.07%, both at 0.4268, so no
% search improvement can help and only a higher tier can.
S.lmit_min = [2000 2000 8000 8000];   S.lmit_min = S.lmit_min(tier+1);

S.tier = tier;   S.tier_base = base;   S.tier_Dup = Dup;   S.dim = dim;
end
