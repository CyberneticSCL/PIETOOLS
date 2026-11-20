# Per-example settings, Phase B (MMP, 11/19/2026; in progress)

Maintainer, 11/19/2026 01:55: "Its fine in Phase B to see which of the default settings tiers works best
as a starting point, but the real goal is to find example-specific changes to these settings which
reduce complexity while maintaining accuracy. For example, the mixed-domain problems benefit from
cross-terms." And 01:45, on the 21 rows no measured setting certifies: "finding settings which work for
these cases is important."

## B0. Tier-level starting points (done 01:37-02:15)

Tags XSP, XSPd0, XLPd0, XSQ, XSQd0 (124 one-variable rows), X2SQ, X2SQd0, X2SP, X2SPd0, X2LPd0 (the
two-variable rows with a Phase A choice), X2HPd0 (gains 45-47), XVP8 (DDE08), XVH (the rows of the 11/17
tree) of `st_run`; results `results_X*.csv` of the 2026-11-18_studies harness. One-variable passes took
133-235 s each. With them `phaseA.pl` (output prefix `phaseB_choices`) picks, per row, the cheapest
setting reaching the best measured result: stripped 'P' with dD 0 is the start of 72 rows; DDE08 reaches
3.9369019 in 'P' at veryheavy (above 'Q' [8 6 6], 3.9144165), so 'P' falls short nowhere.

## B1. The knob search (`ex_search.m`)

From the start of a row, a greedy descent tries one change of one knob at a time, among the knobs the
executive of the row reads, and keeps it when the program has fewer decision variables and the row keeps
its accuracy:
- stability rows: certified under both gates at the library parameter;
- bisection rows: certified at the best threshold less 1e-3 of it (one solve);
- value rows: certified with a value within 1e-3 of the best (of the start value, when the start target
  is looser).

Moves: 'Q' stability - `sop.deg_Q{k}` entries, `dD`, `psatz{k}`, `cap_inc`, `wmin`, `max_raise_Q`;
'P' forms, H2 and synthesis - storage `len_com`, `len_inc`, `len_sos` entries, `patterns`, the same slack
moves, `dw_feas` or `dw_P`, raises, `lpi.var.deg` for synthesis. In two or more variables, after the
reductions converge: cross terms (`wcross` 1, `len_cross{k}` 2) combined with each reduction that failed
alone. The size gate is set below the current size, so a change that does not shrink the program is
refused after its build, before any solve; a change that leaves the size unchanged marks the knob inert
for the row.

Test (EXT, 01:52, four rows; every kept change certified under both gates):

| row | start | final | kept changes |
|---|---|---|---|
| PDE 1 ('P') | 180 / m 24 | 9 / m 5 | storage len_inc{1} 0, len_com{1} 0, len_sos{1} [1;2;2] |
| PDE 33 ('Q') | 409 / 111 (11/18 file st) | 89 / 57 | psatz{2} 'none' |
| DDE10 ('P', gain 0.23644297) | 1,669 / 178 | 1,329 / 157 | len_com{1} 1, len_sos{1} [2;1;1]; gamma 0.23666776 |
| DDE03 ('P', threshold 1.71785) | 440 / 45 | 440 / 45 | none keeps the threshold |

## Runs

- EXS1: every one-variable row (started 01:53).
- EXU: the uncertified rows (`phaseB2_starts.csv`: alternative starts per row - cross terms, the other
  form, heavier tiers, the coercive H2 forms for PDE 50; m capped at 30,000 so no solve passes the
  20-minute limit while the maintainer is online). Maintainer, 11/19: the tests certify several notions
  (only 'P' with strict negativity is exponential stability), so PDE 69 and 70 (documented not
  exponentially stable) are searched, not excluded; undamped beams and waves are stable in the 'Q'
  sense with epneg = 0, so PDE 15 and 16 (Timoshenko beams) are searched in 'Q' (15 also in the primal
  'Q', its library executive being the dual). Not searched: DDE19 and DDF03 (library block defects).

## Patterns for the default tiers (preliminary, 11/19 03:00; `phaseB_tierhints.pl`)

Maintainer, 11/19/2026: "As the results for example-optimized settings come in, look for patterns which
might be used to tweak the default tiers. Obviously, the storage-slack trade for stability/gain problems is
a limitation for default tiers which are built to work on both. After the results are in, we can look at
this more closely." Before the trade pass and the verification, the tuned settings show:
- dD = 0 (slack lift at its minimum) in 89% of the 1-D and 83% of the 2-D stability rows, 76% of the 1-D
  gains, 86% of the H2 norms, 53% of the synthesis rows; the 11/18 study B2 measured that dD = 0 as a
  default loses six 1-D certificates.
- 1-D stability: built on stripped in 58 of 71 rows, form 'P' in 63; the storage cut further (len_com{1} = 0
  in 51%, len_inc{1} = 0 in 32%).
- 1-D gains: stripped in 15 of 17, 'P' in 14; storage changes rare and small (median start/final variables
  1.10 against 2.06 for stability), so the storage sets the gain bound.
- 2-D stability: light in 12 of 18, 'Q' in 14; cap_inc 'auto' 33%, psatz{2} 'none' 28%, deg_Q{2} 0 17%.
- 2-D gains (3 rows): heavy, 'P'; dw_P 1 and a simpler 2-D storage (uniform patterns, smaller len_sos).
- Synthesis (1-D): var.deg{1} [0 2 2] 53%, dw_P 0 37%.

## Open items for later (maintainer's word)

- PDE 8, `PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das` (Das 2018, Example 2): maintainer, 11/19/2026: "Its
  possible there is a typo in this example. Flag it to investigate later, but leave it alone for now."
  As coded the 10/28 audit found max Re lambda = +2.957 (converged in N) against the documented
  stability for lam < pi^2; no setting certifies it (both forms, every tier). The file is not edited and
  the row is not searched. A lead: the older copy `Stability_Examples/PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das`
  differs from the library file in one code line only: the library file sets `PDE_t.x{1}.term{2}.D=1;`,
  the older copy does not, so the two give one term of the PDE different derivative orders.

## Examples outside the PDE library (maintainer 11/19: "Add all 7"; "Don't add duplicates")

Compared with their library counterparts, comments and blanks removed:
- duplicates, not added: `PIETOOLS_PDE_Example_Lamare` (identical code), `_Saba` (GUI file paths only),
  `_Transport_Eq` (the mirror image of case 1: velocity +1 with its condition at s = 1 against -1 at
  s = 0), `_Heat_Eq_with_ODE_Das` (one line, above; Das is left alone);
- added: `PIETOOLS_PDE_Example_Diagne` (gains k1 = 1, k2 = 2 against -0.5, -4.5 in case 3, and another
  layout of the feedback blocks: a different closed loop), `PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin`
  (Robin conditions; unstable for lam > 1.7071, default lam = 1.705), `2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D`
  (the system of case 44 with the delay states on a second variable: a two-variable mixed-domain PIE,
  where case 44 puts them on the variable of the wave).
- EXS2: the two-variable rows, after the B0 two-variable passes; the starts are regenerated first.

State per row: `exs_state/<key>.mat`; logs `exs_<tag>.csv`; results `exs_<tag>_final.csv` (the st
statements from `lpisettings_sop(tier)`); MATLAB logs in `logs/` here.
