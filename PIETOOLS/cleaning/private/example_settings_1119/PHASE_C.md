# Phase C – the example-specific settings written into the examples (11/19/2026)

Maintainer, 11/19/2026: "all examples should have optimized settings files. There is no default."
Phase C writes the settings that Phase B found and verified (PHASE_B.md) into the example files and
the DDE/NDS library blocks, and wires the library run path so that each example runs its own
settings and names the form of each test. All edits are uncommitted, for review (maintainer: "Yes,
uncommitted").

## 1. What was written

| Target | Count | Tool | Content |
|---|---|---|---|
| PDE example files without an `st` | 46 (X1, PDE 53 and 54 after 04:20) | `phaseC_write.pl` | output `st`; OUTPUT help line; note in the description (form, settings, measured result); header log entry with the old signature as a (was) line; the form of each test (`stability_form`, `Hinf_gain_form`, `H2_norm_form`) next to its flag; the `st` block before the GUI branch |
| PDE example files with an 11/17–18 `st` | 11 | `phaseC_update.pl` | changed statements kept as (was) lines; new statements added; dropped statements return to the tier value; note, log entry, form line |
| Indexed files (Lamare 2.x, Saba 4.x, Euler–Bernoulli 13.x) | 3 | by hand | one `st` per case; Euler–Bernoulli sets the form per case (13.1 `'P'`, 13.2 `'Q'`) |
| `PIETOOLS_PDDE_Ex.m` (PDE 50) | 1 | by hand | form `'oco'` (`H2_norm_form`), light with `dD = 0` |
| DDE and NDS/DDF library blocks | 37 (25 DDE, 12 NDS/DDF; 38 edits, DDF01 written twice) | `phaseC_blocks.pl` | `st` lines and forms in each block, at its comment level (the block active in the DDE library, Example 2b, uncommented by hand); 11/18 `st` lines kept as (was) with two comment levels |
| `examples_PDE_library_PIETOOLS.m` | 60 call sites (PDE 53 and 54 added after D145) | `phaseC_library.pl` | `[PDE_t,st] = F(...)`; the run path passes `st` as `settings`; `process_inputs` resets the forms; the 2-D message no longer advises 'stripped' |
| `PIETOOLS_auto_execute_sop.m` | – | by hand | `Hinf_gain_form` and `H2_norm_form` (dual: first letter o → c) |
| `PIETOOLS_DDE.m`, `PIETOOLS_DDF.m` | – | by hand | `Hinf_gain_form`; `PIETOOLS_DDF.m` takes the block `st` as `settings`, as `PIETOOLS_DDE.m` did since 11/18 |

61 example files, 3 libraries and 3 scripts: 67 files, 1,839 insertions, 142 deletions (git diff
--shortstat, 05:50). The edit logs are `phaseC_write_applied.txt`, `phaseC_update_applied.txt`,
`phaseC_blocks_applied.txt` and `phaseC_library_applied.txt`.

## 2. Verification

Two paths, at different parameters:

1. The harness (`st_run` EXV1B, EXV2, EXVX) runs each row of the plan with its planned `st` and
   forms: thresholds by bisection, values for the gains, the dual rows at the catalog parameter.
   EXV1B 115/115 one-variable rows, EXV2 23/25 two-variable rows (the two PDE 50 rows; see §4),
   EXVX 2/2 (X1).
2. The library path (`libcheck.m`, tag LIBCHECK) calls each example file as the library calls it,
   at the defaults of the file, and runs `PIETOOLS_auto_execute_sop` with the `st` and the forms
   the file returns. This path tests the wiring and the defaults the user meets.

Library path (`libcheck_results.csv`; the last line of each example and test counts): 84 tests
over 71 cases in 61 example files (every PDE example file with an `st`, the indexed cases
separately, and X1; PDE 54 added at 05:45), 83 optimal, 1
infeasible: the PDE 50 dual in 'cco' (1.66e6 variables), as its note states. Pass 1 (04:23–04:31)
failed the duals of PDE 7 and 11 and errored on PDE 50; the causes were defects 1 and 2 of §5, and
pass 2 (04:30–04:44) passes PDE 6, 7, 11 and the PDE 50 primal (0.87797144 at 503,633 variables,
577 s with the dual) after the fixes. The executives report on their own flags here, not on the
harness gate |feasratio − 1| <= 0.3; the harness verification above applies that gate to the same
settings.

DDE and NDS/DDF blocks: `blockcheck.m` (tag BLOCKCHECK) uncomments each block with its `st` and
forms in a copy of the library, runs `PIETOOLS_DDE.m` or `PIETOOLS_DDF.m`, and reads back the
settings the executive received, the form, and the verdict by the harness gate. Each library ships
with one block active (DDE: an H2 example, NDS/DDF: DDF03), which the copy comments out, as a user
switching examples does. Result: 60 tests over 35 blocks certified, every one with the settings of
its block (`settings_match` true) and the planned form; the two H2 blocks (DDE25, DDE26), rerun with
the H2 flags (BLOCKCHECK2), certified primal and dual (4 tests). All 37 planned blocks pass.

## 3. Examples added from outside the library (11/19)

| File | Verdict |
|---|---|
| `PIETOOLS_PDE_Ex_Reaction_Diffusion_Robin` (X1) | certified up to lam = 1.70543 (threshold 1.7071) at 189 variables, 'P' stripped with `len_com{1} = 1`, `len_sos{1} = [0;1;1]`; dual certified at 189 (stripped as is: 280) |
| `PIETOOLS_PDE_Example_Diagne` (X2) | infeasible at every tier in 'P' and 'Q'. The file is the predecessor of the library file `PIETOOLS_PDE_Ex_Diagne.m`: its BCs give x_1 (speed lm1 = −3.76) at s = 0 and x_2, x_3 (speeds 0.506, 5.255) at s = 1, the outflow boundaries, so the PDE is not well posed forward in time. The library file reorders the BCs (permutation `Ppm`) and changes the gains to k1 = −0.5, k2 = −4.5. Not edited. |
| `PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D` (X3) | infeasible at stripped, light and heavy, 'P' and 'Q', at mu = 0.375 (documented stable for mu < 0.5); the one-variable form of the same system (PDE 44) fails likewise (§4) |

## 4. Examples without a certified setting

Each row ran the tier passes (stripped to veryheavy, 'P' and 'Q', with and without dD = 0) and
the alternative starts of the search (EXU: cross terms, the other form, slack-focused starts; 3 to
14 per row). The solver signature separates three cases: a certificate of infeasibility
(feasibility ratio −1), an inconclusive exit (MOSEK optimal or unknown at a feasibility ratio far
from ±1), and a program above the harness caps.

| Example | Signature | Diagnosis |
|---|---|---|
| PDE 4.5 (Saba) | infeasible, every tier and form, until 11/19 18:10 | unstable as coded (§7): sig2 = 2.048 for the 1.048 of its source (§10). Corrected 11/19 (maintainer: "Apply all four"); certified at 720 variables (stripped, dD 0), library path optimal |
| PDE 8 (Das) | not searched until 11/19 18:45 | as coded max Re lambda = +2.957 against the documented lam < pi^2: stable only for lam < −5.7799 with the flux input +Bxr of the paper. The sign reversed (maintainer: "probably a typo in the paper. so reverse the sign"): threshold 9.7946; certified up to lam = 9.7925 in 'Q' at 2,496 variables (stripped, deg_Q{1} [2 1 1]), the dual 'Qd' at the default lam at 2,180; library path optimal for both |
| PDE 15, 16 (Timoshenko 2, 3) | infeasible with certificate on 12 of 14 and 7 of 7 starts ('Q', 'P', duals; dD 2, Dmax, wmin [4 4 0], psatz none, slack raises); two 'P' starts of PDE 15 inconclusive (ratio −0.55, −0.31) | unstable as coded, so no certificate exists (§9, 11/19 16:40). PDE 15 drops the bending term E·II·φ_ss of its documented PDE from both input formats (A2 = 0; x3_t = x2 − x4 − c x3), and its header says "Unstable."; at the default c = 0, x(t) = (w1, tφ1, φ1, tφ1) with φ1 = 2s − s², w1 = s² − s³/3 solves the PDE and its BCs (Wolfram), and ‖x4(t)‖² = 8t²/15. PDE 16 (header: "stable"): the characteristic determinant of w_tttt + w_tt − 2w_ttss + w_ssss = 0 under its BCs (w_ss(1) = w(1), w_sss(1) = w_s(1)) has the roots 1.5207 + 1.0639i, 2.2654 + 4.3009i, 2.7484 + 7.5350i (40 digits; Chebyshev collocation at 40 and 56 nodes agrees), so it grows exponentially. Files of other developers (Sachin 2022, Declan 2024, Dan Braghini 05/2026): reported, not edited. PDE 15 corrected 11/19 18:10 (maintainer: "Apply all four"; the bending term of its source, §10): certified in its dual test 'Qd' at veryheavy with dD 0, 74,593 variables ('Qd' at stripped, light, heavy infeasible; the primal 'Q' certifies at light, 21,413, at c = 0, 0.1 and 1); library path optimal. PDE 16 removed from the stability examples 11/19 18:40 (maintainer: "remove it from the stability examples"): the file, its GUI file and its command-line copies with the same BCs (`stability_16.m`, `stability_16_2.m`) in `orphaned_file/private/`; the library entry stops with an error naming the removal |
| PDE 36 (2-D telegraph) | inconclusive (ratio 0.375 and 2.39 at 52,880 and 52,964 variables, MOSEK optimal; unknown at 744,164 and 1,775,200); infeasible at 'P' light (754,180) | exponentially stable (damping 2 lam = 2, Dirichlet BCs); the inconclusive exits point to a weakly feasible program |
| PDE 39 (2-D Navier–Stokes, Antonelli) | inconclusive (ratio 2.64 at 316,737; uncapped, ratio 7.14 at 648,900, m 30,721, MOSEK optimal, eta 2.4e-9, 271 s); infeasible ('P', 238,905) | open |
| PDE 41 (heat, delayed boundary input) | inconclusive on all 12 starts: MOSEK unknown on 11 and optimal on 1, feasibility ratio between −0.80 and 0.43 | open; heat with lam = 1 < pi^2 and a transport state with zero inflow when u = 0 |
| PDE 44 and X3 (wave, partly delayed boundary damping, second-order form) | infeasible (feasibility ratio about −1) on 8 of 9 starts (27,731 to 1.56e6 variables) | exponentially stable at the tested mu: spectral abscissa −0.32 (mu = 0.375) and −0.29 (default mu = 0.4) from cosh lam + k sinh lam ((1−mu) + mu e^{−lam tau}) = 0, k = 1, tau = 1 (Wolfram, 11/19). The same system in Riemann form (PDE 59, states w_s, w_t and the delayed w_t) certifies at 136 variables. The form of 44 (states u, u_t, an ODE state x1 tied to the BCs, and u(1, t − tau s) beside u_t(1, t − tau s)) defeats the LPI, not the system; a reformulation as in PDE 59 is an example-design decision |
| PDE 50, dual | 'cco' infeasible (96,417 and 1.66e6 variables); 'c' unknown at every tier (values 305 to 1.3e5) | the primal certifies in 'oco' (0.877971); the dual stays open |
| PDE 56 (2-D heat, 3-D delay state) | uncapped light 'Q' with dD 0: infeasible (840,156 variables, m 33,352, 178 s) | the D2 program of the 11/18 studies (D119, D120: capped light infeasible; the uncapped program a Sol job held for the maintainer, outside the Sol approval for PDE 53–55) |
| PDE 55 (3-D coupled heat, lam = 1 against 16.2) | 'Q' light or stripped (dD 0, deg_Q{3} 0): 12,959,852 variables, m 189,509; with cap_inc 'auto' the same (dD 0) or 29,455,660, m 277,565 (dD 1); 'P' stripped: dD 0 infeasible (114,704), dD 0 with wmin [0 0 1] the same program, cap_inc 'auto' infeasible (1,808), dD 1 7,340,048 variables, m 105,717; 'P' light 7.36e6 (dD 0) and 2.80e7 | held: 'P' stripped with dD 1 certifies the one-field PDE 53 (580,612 variables, m 11,753); on PDE 55 it needs a Sol node, projected from the PDE 54 solve (m 32,846, 55 s per iteration on 7 threads) by m^3 to about 8 min per iteration on 64 threads, 3-4 h, 100-200 GB, about 400 CHE, above the pre-flights approved so far (10/30: 45-90 CHE). Approved 11/19 16:20; Sol job 2067511 (public, 64 cores, 256 GB, 6 h limit, at most 768 CHE; the working tree of 16:28, one build and one solve through ex_search, `sol55.m`), queued at 16:30 with a start estimated for 11/20 05:43 |
| PDE 65 (Turing pair at 0.9 d*) | infeasible (499,640 and 1.75e6; uncapped light 'Q' with dD 0, deg_Q{2} [1 0 0], wcross 1: infeasible at 2,039,468 variables, m 28,303, 673 s) | open; 10% below the Turing threshold |
| PDE 69, 70 (wave–heat 2-D, Beale acoustic) | infeasible with certificate on every start (336,681 to 1.01e6 variables) | documented not exponentially (69) and not uniformly (70) stable; with epneg = 0 the certificate needs a non-increasing storage |
| DDE19, DDF03 | inaccurate or unknown at every tier until 11/19 17:30 | no regulated output reached the converter; with their outputs (§7) both certify: DDE19 1.963986 at 265 variables, DDF03 0.439636 at 4,673, each through PIETOOLS_DDE.m / PIETOOLS_DDF.m with its st |
| X2 | infeasible, every tier and form | ill-posed predecessor of the library Diagne file (§3) |

The signatures do not depend on the solver: with SeDuMi in place of MOSEK (EXSED, 05:50), PDE 15, 16
and 44 are infeasible at feasibility ratios −0.54 to −1.18 ('Q' light and 'P' stripped), PDE 41
inaccurate at 0.32 and 0.37, and PDE 36 unknown at 0.64.

Two measured pairs bear on the formulation of these examples. PDE 44 and X3 (wave with partly
delayed boundary damping, states u, u_t) fail on every start, while PDE 59 (the same system in
Riemann form, states w_s, w_t) certifies at 136 variables. PDE 15 (Timoshenko, states w_t, w_s,
phi_t, phi, damping c = 0 by default) fails on every start, while PDE 14 (Timoshenko, states w_t,
the shear force, phi_t, the bending moment) certifies at 144 variables. A position state is not
fatal in general: PDE 38 (2-D wave, states x, x_t) certifies at 380,804 variables. Rewriting an
example in energy states is an example-design decision, left to the maintainer.

## 5. Defects found and fixed in Phase C

1. `phaseC_plan.pl` chose the stability form of an example by a vote over its rows and broke a
   primal/dual tie by Perl hash order, independently of the row whose `st` it used. Two runs could
   name different forms; the files of PDE 6, 7 and 11 received the form EXV1B had not verified.
   Fixed: the form of the largest row of each kind; the three files now carry the EXV1B forms
   (6 'P', 7 'Q', 11 'P'); DDF01 keeps the verified 'P' through an override.
2. `PIETOOLS_auto_execute_sop.m` (11/19 edit) derived the H2 dual form by `strrep(hform,'o','c')`,
   which maps 'oco' to 'ccc'; the library check of PDE 50 errored. Fixed: `['c',hform(2:end)]`.
3. Hand edits through the Edit tool dropped trailing blanks of four unchanged lines; restored, so
   the diff holds no whitespace-only changes.
4. `phaseC_blocks.pl` wrote every block one comment level down; the block active in the DDE
   library (Example 2b, H2 norm) is uncommented, so its `st` lines stayed comments when the
   shipped `PIETOOLS_DDE.m` runs it. Set to level 0 by hand (05:00); the block check passes it.
5. A first version of `blockcheck.m` left the active block of each library running after the
   block under test, so `prog` came from the H2 example; fixed (the copy comments it out) and rerun.
6. `PIETOOLS_DDF.m` (11/12 edit) set `settings.lpi.pos.sep = 1` before the H-infinity controller and
   estimator syntheses; since 11/17 every tier writes `lpi.pos` as a storage row (`len_com`), which
   `exec_sop_settings` refuses with stock fields, so both syntheses stopped with an error for any
   settings. `PIETOOLS_DDE.m` sets no `sep`. Fixed 11/19 17:38: `sep` only for a stock-format
   `lpi.pos`. `blockcheck.m` also matched only commented blocks; it now tries the lines as they are
   for the active block (DDF03).

## 6. Patterns for the default tiers (for review; no tier file changed)

Measured on the final settings of the search (`phaseB_tierhints.pl`, `phaseB_patterns.pl`) and on
the paired tier passes with and without the slack lift degree (stripped tier, every row; results
files `XSP`/`XSPd0`, `XSQ`/`XSQd0`, `X2SP`/`X2SPd0`, `X2SQ`/`X2SQd0`).

1. **The slack lift degree dD.** The final settings carry `dD = 0` in 86% of the one-variable
   stability rows (73), 80% of the two-variable stability rows (20), 76% of the one-variable gain
   rows (17), 86% of the H2 rows (7) and 53% of the synthesis rows (19). In the paired passes:

   | Pass pair | Rows | Pass with both | Result equal | Pass with dD = 1 only | dD = 0 only | Variables dD 1 / dD 0 (geometric mean, median, both passing) | dD 0, then dD 1 on failure / dD 1 always (summed variables) |
   |---|---|---|---|---|---|---|---|
   | 1-D 'P' | 124 | 87 | 80 | 17 | 2 | 1.60, 1.56 | 0.77 |
   | 1-D 'Q' | 124 | 90 | 84 | 10 | 5 | 1.37, 1.42 | 1.02 |
   | 2-D 'P' | 14 | 7 | 7 | 4 | 1 | 2.73, 2.00 | 1.10 |
   | 2-D 'Q' | 12 | 9 | 9 | 1 | 1 | 2.54, 3.09 | 0.98 |

   Rows passing with dD = 0 only: DDE26 dual and PDE 20 dual ('P'); DDE12, DDE26 dual, NDS11,
   PDE 22, 68 ('Q', one variable); PDE 42 ('P', two); PDE 38 ('Q', two). Where both pass, dD = 0
   costs 1.4–2.7× fewer variables at an equal result on 80 of 87, 84 of 90, 7 of 7 and 9 of 9
   rows. A tier that tries dD = 0 first and raises dD on failure, with the failed attempt counted,
   sums to 0.77× the variables of dD = 1 in one-variable 'P' and 0.98–1.10× elsewhere: the gain
   holds for one-variable 'P' only, unless the failure can be predicted. The executives raise
   degrees, not dD, at present. (A list-context error in the first count of 04:50 reported no row
   passing with dD = 0 only; corrected 05:40.)
2. **Form.** One-variable stability: 'P' on 65 of 73 rows; one-variable gains: 'P' on 14 of 17;
   two-variable stability: 'Q' on 16 of 20; two-variable gains: 'P' on 4 of 4.
3. **Tier.** One-variable rows settle on stripped (stability 60 of 73, gains 15 of 17);
   two-variable stability on light (13 of 20).
4. **Storage against slack** (kept moves / moves tried, variables removed). Storage cuts keep the
   verdict on 251 of 471 stability moves (78,213 variables) and the accuracy on 32 of 198 gain moves
   (114,265); slack cuts keep 25 of 76 stability moves (300,620) and 3 of 48 gain moves (1,195,110).
   The slack holds most of the variables; storage cuts are frequent and small, and on gains they
   lose accuracy in 84% of the trials. Trades (one family raised to cut the other) are kept 18 times
   in 231 trials: gains 2 of 59, stability 5 of 71, synthesis 11 of 101. Against the hypothesis of
   11/19 02:45 (storage matters for gains, slack for stability): on stability rows slack cuts lose
   the verdict more often than storage cuts (51 of 76, 67%, against 220 of 471, 47%), as it
   predicts; on gain rows both families lose accuracy in most trials (storage 166 of 198, 84%;
   slack 45 of 48, 94%), so these data do not separate them.
5. **Two-variable specifics.** `deg_Q{2}` lowered to `[0 1 1]` or `0` on 35% of the two-variable
   stability rows, `psatz{2} = 'none'` on 25%, `cap_inc = 'auto'` on 30%; the two-variable gains use
   `patterns{2} = 'uniform'` on 3 of 4.
6. **Synthesis.** `var.deg{1} = [0 2 2]` on 53% of the synthesis rows.
7. **`deg_Q{2} = [0 1 1]` as a two-variable default, measured** (TQ2ctl against TQ2011, 06:00:
   light 'Q' on the 19 two-variable stability rows, as is and with that degree). Light certifies 16
   rows (38 unknown, 58 infeasible, 61 above the 2M size gate at 2,430,327 variables); with
   `[0 1 1]` 11 of them stay certified and 5 are lost (PDE 43, 57, 60, 62, 63). Variables drop on
   PDE 33 (52,660 → 18,081), 42 (69,915 → 42,261) and 59 (5,214 → 3,874), stay within 0.1% on 34,
   35, 37, 40, 64, 67, 68, and grow on 66 (56,548 → 187,281). The change pays on a few examples and
   is no candidate for the tier.
8. **Other frequent changes as defaults, measured the same way** (`pair_compare.pl`, 06:00–06:30):

   | Change (pass) | Rows | Pass both, result equal | Lost | Gained | Variables control/test, geometric mean |
   |---|---|---|---|---|---|
   | `cap_inc = 'auto'`, two variables, light 'Q' (TQ2cap) | 19 | 15, 15 | PDE 60 | PDE 61 (light: 2,430,327, above the size gate) | 1.54 |
   | `psatz{2} = 'none'`, same rows (TQ2psn) | 19 | 6, 6 | 10 (PDE 34, 35, 37, 42, 43, 57, 62, 63, 64, 66) | PDE 38 | 1.45 |
   | `len_com{1} = 0`, one variable, stripped 'P' (T1Pc0 against XSP) | 124 | 95, 82 | 9 | 2 | 1.17 |
   | `len_inc{1} = 0`, same (T1Pi0) | 124 | 85, 69 | 19 | 1 | 1.17 |

   The triangular cap is the one candidate among them: on the two-variable stability rows at light
   it keeps 15 of 16 certificates with equal results, gains PDE 61, and cuts the variables 1.54×
   (geometric mean); it loses PDE 60. In the storage form 'P' the cap collapses the slack instead:
   on the two-variable gain rows at light 'P' (TG2cap against TG2ctl) the programs shrink from
   200,986–203,522 to 1,851–2,226 variables and all four become infeasible (PDE 45, 46, 47, 71), as
   PDE 55 did in 'P' stripped (1,808). The candidate is the direct form 'Q' only. At heavy 'Q'
   (TQ2hcap against TQ2hctl, 06:30–06:48) the cap keeps all 19 two-variable stability certificates
   and cuts the variables on 9 of them (PDE 42 249,298 → 146,674, 43 215,188 → 126,228, 57 248,497
   → 151,473, 58 302,957 → 182,429, 59 10,581 → 6,225, 60 4,700 → 2,764, 61 841,042 → 575,082, 67
   34,900 → 32,300, 68 12,564 → 4,392), equal on the other 10: 0.89× the summed variables.
   Candidate for the maintainer: `cap_inc = 'auto'` in the two-variable entry of the direct form at
   light and heavy (light: 15 of 16 kept, PDE 60 lost, PDE 61 gained; heavy: 19 of 19). On the one-variable rows at light it was measured on 11/18
   (D102): −17% variables, 6 rows better and 8 worse, and kept off in every tier; the two-variable
   entry was not measured then. `psatz{2} = 'none'`, `len_com{1} = 0` and `len_inc{1} = 0` lose
   certificates or accuracy on many rows (13 and 16 changed results on the rows the two storage
   cuts pass, e.g. DDE01 2.04 → 1.81 and 2.04 → 2.006).

## 7. Findings for the maintainer

- Saba 4.5 as coded (sig2 = 2.048) is unstable: phi' = M(lam) phi with
  M = [−r1 lam, r1 sig1; −r2 sig2, r2 lam] and the BCs give F(lam) = [−pb 1] expm(M(lam)) [qb; 1];
  F has the real root lam = 1.682, and F(0) = 0 at sig2 = 1.0484, the "max sig2=1.049" of the comment
  of the file (Wolfram, 11/19/2026). At sig2 = 1.048 the largest real part is −0.0017. A typo
  (2.048 for 1.048) would explain it; the parameter is left as it is, the note of the file states it.
- PDE 8 (Das): left alone, as instructed.
- DDE19 (wind tunnel, `Hinf_control`) names its regulated output with the fields `C0`, `C{1}`,
  `D1`, `D2`; `initialize_PIETOOLS_DDE` dumps the fields to variables and reads `C1`, `C1i`,
  `C1di`, `D11`, `D12`, ..., zero-filling those absent, so z is identically zero. The gain bound
  then tends to 0: every tier returns 5e-4 to 1.4e-3 with MOSEK inaccurate or unknown, against the
  documented 1.964 (degrees 2,4 and 4,8) and 1.9369 (6,12). No `st` written; the block needs its
  field names updated (not edited here).
  **Resolved 11/19 17:40** (maintainer: "find out if they have an established regulated output.
  Otherwise ... add one and document it"). Neither the model source [5] (Manitius 1984, pole
  placement) nor the arXiv version of [4] (1806.08071, whose examples do not include the wind
  tunnel) defines z for this model. The block has carried z = [x1; x2; .1u] since PIETOOLS 2021b,
  the form of the examples of [4] (z = [x; .1u]); it now sets it with `C1`, `C1i`, `D11`, `D12`,
  documented in the block. Certified at every tier: gamma 1.963985 (stripped 1,252 variables, light
  1,876, heavy 2,673), the documented 1.964; the documented 1.9369 (degrees 6,12) is not reproduced
  by any tier. The search EXSCTL cut stripped to 265 variables at the same bound (`st` in the
  block: len_inc 0, len_com 0, len_sos [1;2;2], dD 0, var.deg [0 2 2]); PIETOOLS_DDE.m certifies it
  (1.9639856, settings matching). The model of [5] has a discrete delay (kappa m' + m =
  k theta(t − tau0), eq. (6.7) of Mazanti, Boussaada, Niculescu, arXiv 2002.06128); the block
  declares `DDE.Adi{1}`, a distributed delay, since 63831988 (12/2022; `Ad{1}` before, outside the
  struct). With `DDE.Ai{1}` the bound is 1.963867 at every tier. Reported in a note of the block, not
  changed. **Applied 11/19 18:05** (maintainer: "Apply all four"): `DDE.Ai{1}`. The st of the
  distributed delay (265 variables) ends unknown there (ratio 0.594); the search EXSFIX2 gives 221
  variables (stripped; len_inc 0, len_com 0, len_sos [1;2;2], dD 0, var.deg [0 2 2], dw_P 0), gamma
  1.963870; PIETOOLS_DDE.m certifies it (1.9638696, ratio 0.751, settings matching).
- DDF03 (showering people, `Hinf_control`) leaves all three output options commented out and sets
  only `Bv`, `Cr`, `Cv`: no regulated output, the same collapse (5e-4 to 1.4e-3, inaccurate or
  unknown). No `st` written; not edited.
  **Resolved 11/19 17:40.** The source, Peet, Automatica 2021 (arXiv 1910.03881), eq. (15), sets
  z = [sum_i T1i; .1 sum_i ui], option 2 of the block, now active and documented. Certified at every
  tier: gamma 0.447532 (stripped, 3,557 variables), 0.442635 (light, 5,669), 0.440633 (heavy,
  8,465). The search EXSCTL found 0.439636 at 4,673 variables (`st` in the block: heavy with len_sos
  [4;2;2], cap_inc 'auto', var.deg [1 2 2]). The source sets tau_i = i; the block sets tau_i = n
  (its comment says i); with tau_i = i the bounds are 0.425834, 0.423643, 0.421337. Reported in a
  note of the block, not changed. Through PIETOOLS_DDF.m the block first stopped in
  exec_sop_settings: the controller and estimator branches set `lpi.pos.sep = 1`, which a storage
  row rejects, and every tier is a storage row since 11/17, so both syntheses failed for any
  settings. Fixed (§5, item 6); PIETOOLS_DDF.m then certifies the block (0.43963563, settings
  matching). **Applied 11/19 18:05** (maintainer: "Apply all four"): tau_i = i. The search EXSFIX keeps
  the st (4,673 variables) at gamma 0.42054329 (heavy as it is: 0.42133672 at 8,465); PIETOOLS_DDF.m
  certifies it (settings matching). DDE18, the DDE form of the same network, also sets tau_i = n
  (not changed).

## 8. Second round (11/19 04:50–06:30)

- **Resumed searches** (EXS2X, EXS2Y): the ten two-variable rows whose first search ended on the
  40-build budget, 30 more builds each from the saved state. PDE 57, 58, 61, 63, 64 converged and
  PDE 45, 46, 47, 62 reached the extended budget without a further cut. PDE 71 (H-infinity gain,
  exact gamma* = 0.04571818) dropped from 202,922 to 94,987 variables by the cross weight
  (`wcross = 1`) with `len_com{2} = [2 1]`: gamma 0.0459804 (+0.574 %) against 0.0459442 (+0.494 %),
  within the search tolerance (1e-3 of the best measured). The file carries the new st; the
  library path gives 0.045980398 at 94,987 (optimal).
- **3-D.** 'P' stripped as it is certifies PDE 53 at 580,612 variables, m 11,753 (112 s), against
  668,443, m 18,102 for the light 'Q' setting written at 04:20, and PDE 54 at 1,835,012, m 26,429
  (671 s, feasibility ratio 0.957), where light 'Q' (1,922,843, m 32,846) ended optimal at ratio
  1.69, outside the gate. The searches from it (EXS53, EXS54) find no cut: dD 0 is infeasible
  (28,676) on both, every other knob inert. Both files now return 'P' stripped; PDE 54 joins the
  certified examples. PDE 55 (two fields) needs 7,340,048 variables, m 105,717 in this setting: held
  for a Sol decision (§4).

## 9. Diagnostics of the uncertified rows (11/19 16:30–17:00)

`diag_amt.m` (st_run tags DIAGA1–DIAGA7) builds the stability LPI of a row with the generator A
replaced, at the settings of a failing start, and solves it:

- **A0** (A = 0): only the positivity constraints of the form act. Infeasible means the storage
  parameterization admits no coercive storage, whatever the dynamics.
- **AmT** (A = −T, so x_t = −x): infeasible means the negativity slack cannot carry 2T'T ('Q') or
  2T'PT ('P').
- **A − cT**: certified exactly when a storage of the program gives dV/dt ≤ 2cV under A, a growth
  bound c. A − T can cancel terms of A; where the program at c = 1 is smaller than the program of A
  (PDE 15, 44), it is not the same program.

Light tier; decision variables in parentheses; the A column of the two-variable rows is from the
EXU logs.

| Row | Form, start | A | A0 | AmT | A − cT |
|---|---|---|---|---|---|
| 15 | 'Q' | infeasible (73,364) | certified (5,520) | certified (10,608) | c = 1 certified (20,516) |
| 15 | 'P' | infeasible (16,000) | certified (1,680) | certified (4,480) | c = 1 certified (4,480) |
| 16 | 'Q' | infeasible (119,514) | certified (22,253) | certified (44,074) | c = 1 infeasible (145,098) |
| 16 | 'P' | infeasible (17,225) | certified (1,680) | infeasible (17,225); certified with dD 2 (28,825) or Dmax (23,453) | c = 1 infeasible |
| 41 | 'Q' | unknown, ratio 0.08 (19,776) | certified (1,380) | certified (2,652) | c = 1 unknown, ratio −0.68 |
| 41 | 'P' | unknown, ratio −0.10 (4,000) | certified (420) | certified (1,120) | c = 1 unknown, ratio 0.01 |
| 44 | 'Q' | infeasible (93,828) | certified (14,093) | certified (27,677) | c = 1 certified (30,320); c = 0.5, 0.25, 0.1, 0.05 infeasible (93,828); heavy builds a program of the same size (93,828, m 3,648): A, c = 0.25, 0.1 infeasible |
| 44 | 'P' | infeasible (16,419) | certified (1,787) | certified (4,707) | c = 1 infeasible |
| X3 | 'Q' | infeasible (28,150) | certified (4,123) | certified (7,867) | c = 1, 2 infeasible |
| X3 | 'P' | infeasible (32,406) | certified (2,174) | certified (16,614) | c = 1 infeasible |
| 69 | 'Q'; dD 0, deg_Q{2} 0, wcross 1 | infeasible (426,641) | infeasible (334,341) | infeasible (380,241) | – |
| 69 | 'Q'; dD 0, cap_inc 'auto', wcross 1 | infeasible (1,005,460) | certified (47,025) | certified (93,060) | c = 1 infeasible (1,005,460) |
| 65 | 'Q'; dD 0, deg_Q{2} 0, wcross 1 | infeasible (499,640) | infeasible (479,992) | infeasible (479,992) | – |
| 70 | 'Q'; dD 0, wcross 1 | infeasible (494,572) | certified (32,700) | certified (64,600) | c = 1, 2, 4 infeasible (494,572) |
| 56 | 'Q'; dD 0 | infeasible (840,156) | certified (259,576) | certified (517,596) | – |
| 36 | 'Q'; dD 0, wcross 1 | unknown, ratio 0.375 (52,880) | certified (20,900) | certified (41,360) | c = 1 unknown, ratio 3.23 |
| 39 | 'Q'; dD 0, deg_Q{2} [2 1 1], wcross 1 | unknown, ratio 2.64 (316,737) | certified (47,214) | certified (93,249) | c = 1 unknown, ratio 3.17 |

Readings:

1. PDE 15 and 16 are unstable as coded (§4). The program agrees: it certifies PDE 15 shifted by
   rate 1 (the growth is linear) and not PDE 16 (the spectrum reaches real part 2.75).
2. The EXU starts of PDE 69 and 65 with deg_Q{2} = 0 admit no coercive storage (A0 infeasible), so
   their infeasibility says nothing about the dynamics. At the default deg_Q, PDE 69 passes A0 and
   AmT. A0 is a check of the storage side of a 'Q' start, independent of the row.
3. For every other row and form, the storage side and the negativity slack carry trivial dynamics;
   the obstacle lies in the dynamics of the row. The exception is PDE 16 in 'P', whose light slack
   lift cannot carry 2T'PT for its fourth-order T; dD 2 or Dmax can.
4. PDE 44: no shift from 0.5 to 0.05 certifies in light or heavy 'Q', against a spectral abscissa of
   −0.29; heavy builds a program of the same size as light, so the tier difference does not reach
   this row. PDE 59 (Riemann form) certifies the same system at 136 variables.
5. PDE 41, 36 and 39 end unknown with or without the shift, while A0 and AmT certify: the
   difficulty is numerical and specific to their dynamics, not a missing decay margin. PDE 41 stays
   unknown at lam = 0 (ratio 0.39) and at tau = 0.5 and 2 (−0.15, 0.34), the same program (19,776),
   so neither the reaction nor the delay length causes it. What remains is the coupling of the heat
   state to the delayed boundary input (inferred by elimination, not isolated).
6. PDE 69 (default deg_Q) and 70 keep the program of A under the shift and stay infeasible at
   c = 1 (69) and at c = 1, 2, 4 (70): at these settings the program cannot bound their growth by
   rate 4, while both pass A0 and AmT. Since T is compact, a large shift need not make the shifted
   LPI feasible for a fixed storage parameterization, so this measures the parameterization, not the
   rows.
7. PDE 56 at dD 0 passes A0 and AmT; its dynamics fail there. The held Sol job (light 'Q' at the
   default dD, a larger slack) is not ruled out by this check.
8. Cross terms on the delay rows (DIAGA9: `wcross = 1`, and `len_cross{2} = 2` with it; light and
   heavy; 'Q' and 'P'), never tried on them before: PDE 41 and 44 build the same programs as without
   them (the knobs do not act on these PIEs) and end as before; X3 grows to 58,915–294,603
   variables and stays infeasible.

## 10. Unstable-as-coded examples: literature check (11/19)

Maintainer, 11/19/2026: "task a subagent with looking at the unstable as coded examples. Looking at
the literature to see if they are not formulated correctly and, if they are unjustified, removing
them". Classes: faithful to a source that claims stability (contradiction reported, no change);
mis-transcription of a justified system (correction reported, no edit – files of other developers);
unjustified, including a superseded duplicate (moved to `cleaning/private/dead_code/`). Sources: arXiv
abstract and ar5iv pages, search snippets, and the local copies in
`GitHub/CyberneticSCL.github.io/Publications/` (`2019/Peet_CDPS_2019.pdf`, `2020/Peet_AUT_2020b.pdf`,
`2024/shivakumar_TAC_2024a.pdf`); no downloads. Scripts: `litcheck_1119/`. Authors from
`git log --follow`.

| Example | Source | Verdict | Action |
|---|---|---|---|
| 15 | [8] §8.1.0.2, mixed hyperbolic-diffusive form | mis-transcription: bending term dropped | corrected 11/19 18:05 (maintainer: "Apply all four"), certified ('Qd', 74,593) |
| 16 | [8] §8.1.0.2 beam as one equation; BCs as in TAC 2025 §VII-B | mis-transcription: BCs | removed from the stability examples 11/19 18:40 (maintainer); file in `orphaned_file/private/` |
| 4.5 | Peet 2021 §10 Ex. 5 (Saba 2020 Fig. 5) | mis-transcription: sig2 2.048 for 1.048 | corrected 11/19 18:05 (maintainer: "Apply all four"), certified ('P', 720) |
| 8 | Das 2018 §6.2 | faithful to the data; claim contradicted | flux sign reversed 11/19 18:45 (maintainer: "probably a typo in the paper"), certified ('Q', 'Qd') |
| X2 | Diagne 2012; superseded by case 3 | unjustified: duplicate, ill-posed | moved |

**PDE 15**, `PIETOOLS_PDE_Ex_Timoshenko_Beam_2.m` (Declan 2022-08-20 e21afbef; Sachin 2022-11-22,
2022-12-04, 2022-12-07; Dan Braghini 2022-12-03, 2026-05-31; Declan 2024-12).
- Source: Peet, Shivakumar, Das, Weiland, "Discussion paper: A new mathematical framework for
  representation and analysis of coupled PDEs", IFAC-PapersOnLine 52(2):132–137, 2019 (library [8],
  "Example 8.1.0.2"), §8.1, paragraph 8.1.0.2. After the energy-state form (case 14) the paper takes
  u2 = w_s, u4 = phi, a mixed hyperbolic-diffusive form with A0 = [0 0 0 0; 0 0 0 0; 0 1 0 −1;
  0 0 1 0], A1 = [0 1 0 −1; 1 0 0 0; 0 0 0 0; 0 0 0 0], A2 = [0; 0; 1; 0] on u4_ss and the five BCs
  of the file. It reports that this form does not appear stable in its states, adds the damping
  −c·u3 to the u3 equation, and reports stability for every c > 0.
- Comparison: A0, A1, B and c agree; A2 does not. The file has `PDE_b.A2 = [0;0;0;0]` and no term
  x = 4, D = 2 in `PDE_t.x{3}`. The 2021b library had A2 = [0;0;1;0] in the batch format; the term
  format (DJ, 05/2021) left the term out, except in its commented-out alternative grouping, and
  72de7c6b (Sachin, 12/07/2022) set the batch A2 to zero. `command_line_input_Examples/stability_15.m`
  lacks the term as well.
- Verdict: mis-transcription (dropped term) of a justified system. "Unstable." matches the use in the
  source: c = 0 is the variant on which the test fails, c > 0 the stable one.
- Correction (not applied): `PDE_b.A2 = [0;0;1;0];`, `PDE_t.x{3}.term{2}.x = 4;`,
  `PDE_t.x{3}.term{2}.D = 2;` (x3_t = x2 − x4 + x4_ss − c·x3), the same term in the header and in
  library line 303; `diff(x3,t)==x2-x4+diff(x4,s,2)` in `stability_15.m`. With the term and c = 0 the
  system is the clamped-free beam: the characteristic function of the two-equation beam has the
  zeros i·omega, omega = 1.27396, 1.89396, 4.64411, 4.83045, 7.81368, 7.92542, 10.96691, 11.04672
  below 12, and none in 0.001 ≤ Re ≤ 6, |Im| ≤ 25 (Wolfram, argument principle, 30 digits): stable,
  not exponentially stable, as Peet, Automatica 125 (2021) 109473, §11.1, states for both beams.

**PDE 16**, `PIETOOLS_PDE_Ex_Timoshenko_Beam_3.m` (Declan 2022-08-20 e21afbef; Sachin 2022-11-22;
Declan 2024-12; Dan Braghini 2026-05-31).
- Sources: the header rewrites the beam of [8] §8.1.0.2 (clamped at s = 0, free at s = 1) as one
  equation; [8] has no one-equation form. Shivakumar, Das, Weiland, Peet, "Extension of the partial
  integral equation representation to GPDE input-output systems", IEEE TAC 70(5):3240–3255, 2025
  (arXiv 2205.03735), §VII-B, has the same equation with rho = A = I = kappa = G = E = 1 and the
  coded BCs w(0) = 0, w_s(0) = 0, w_ss(1) = w(1), w_sss(1) = w_s(1); it simulates the response to
  d(t) = sin t (Fig. 7) and states no stability property.
- Comparison: (alp, bet, gam) = (1, 1, 2) is the beam with unit parameters. Eliminating phi from the
  beam BCs (unit parameters) gives w(0) = 0, w_sss(0) − w_tts(0) + w_s(0) = 0 (phi(0) = 0), w_ss(1) − w_tt(1) = 0
  (phi_s(1) = 0) and alp·w_sss(1) − gam·w_tts(1) = 0 (w_s(1) = phi(1)): the BCs of
  `command_line_input_Examples/stability_16_3.m` (Dan Braghini, 06/19/2026). The coded w_s(0) = 0 is
  the Euler–Bernoulli clamp; the coded conditions at s = 1 have w and w_s where the beam has w_tt and
  w_tts. Measured (characteristic function of w = exp(lam t) W(s), Wolfram, 30 digits): with the
  BCs of `stability_16_3.m` the zeros below |lam| = 12 are the eight i·omega of the two-equation beam
  (above), and none lies in 0.001 ≤ Re ≤ 50, |Im| ≤ 50; with the coded BCs, 8 zeros lie in
  0.001 ≤ Re ≤ 6, |Im| ≤ 12 (1.5207 ± 1.0639i first, §4).
- Verdict: mis-transcription (wrong BCs) of a justified system; "stable" holds for the beam in the
  sense of Peet 2021 §11.1 (stable, not exponentially stable). The coded BCs have a published
  source, TAC 2025 §VII-B, which claims no stability.
- Correction (not applied): in the term format, BC 2 `x4_sss(0) − x3_s(0) + x4_s(0) = 0`, BC 3
  `x4_ss(1) − x3(1) = 0`, BC 4 `alp*x4_sss(1) − gam*x3_s(1) = 0`, and BCs 5–8 (time derivatives of the
  old BC 2) replaced by `x3(0) = 0` and `x4_s(0) + bet*int(x3,s,0,1) = 0`, as in `stability_16_3.m`
  (unit parameters; there x1 = w_t and x2 = w_ttt, here the reverse). The integral BC holds for every
  eigenfunction of the other four (integrate the equation over [0, 1]; lam^2 = −1 is not an
  eigenvalue).

**Saba 4.5**, `PIETOOLS_PDE_Ex_Saba.m` (Declan 2022-08-20 e21afbef, the --follow log reaches the
predecessor 6a5d184e; Sachin 2022-11-22; Declan 2024-12; the 4.5 line unchanged since the 2020a
upload c587fd94, 07/06/2020).
- Sources: Bou Saba, Bribiesca Argomedo, Auriol, Di Loreto, Di Meglio, IEEE TAC 65(7):2941–2956, 2020
  (library [3]), Fig. 5, the stability region (not read: publisher paywall, HAL blocks automated
  access). Peet, Automatica 125 (2021) 109473 (arXiv 1812.06794v4), §10, Example 5: the system of the
  file with r2 = 1.1, sigma1 = 1, q = 1.2 (r1 printed as 0; the file and 4.2–4.4 use 0.5), stability
  verified for the stable region of Saba Fig. 5 by gridding sigma2 and rho; at rho = −0.4, proven for
  sigma2 ≤ 1.048.
- Comparison: 4.5 has pb = rho = −0.4, qb = 1.2, sig1 = 1, r2 = 1.1, and sig2 = 2.048 with the comment
  "max sig2=1.049". Measured (§7): F(0) = 0 at sig2 = 1.0484; largest real part −0.0017 at 1.048;
  real eigenvalue 1.682 at 2.048.
- Verdict: mis-transcription (typo, 2.048 for 1.048). Correction (not applied): `sig2=1.048;` on the
  4.5 line. The unused duplicate `PIETOOLS_PDE_Example_Saba.m` carries the same value.

**PDE 8**, `PIETOOLS_PDE_Ex_Heat_Eq_with_ODE_Das.m` (Declan 2022-08-20 e21afbef; Sachin 2022-11-22,
2022-12-04; Dan Braghini 2022-11-22, 2022-12-01; Declan 2024-12).
- Source: Das, Shivakumar, Weiland, Peet, "Representation and stability analysis of PDE-ODE coupled
  systems", arXiv 1812.07186 (v3, 12/20/2018; library [7], "Example 2"), §6.2: A (4×4), B (4×2), C
  (2×4) as in the file; z_t = lam·z + z_ss + B1·u with B1 = 1; case (a) z(0) = z(1) = 0; the flux at
  s = 0 drives the ODE (v = y) and the ODE output drives the PDE (u = w), L = [0 I; I 0]. The paper
  states that the method verifies stability for lam < pi^2 (case a) and lam ≤ 2.467 (case b,
  z_s(1) = 0), the bounds of the PDE alone, and reports no computed value.
- Comparison: the paper writes z as a scalar, while B and C need two PDE outputs and inputs; the
  file takes z in R^2 with B1 = I and v = z_s(0), the reading consistent with the dimensions.
  Measured (Chebyshev collocation, N = 40, 56, 72 agree to 5 digits; `litcheck_1119/das_*.m`):
  eig(A) = −0.680 ± 5.866i, −1.279, −1.820; as coded, largest real part 0.495 at lam = 0 and 2.957 at
  the default lam = pi^2 − 1, stable only for lam < −5.7799; with the flux sign reversed (−z_s(0), or
  z_s(1)), stable for lam < 9.7946 (−0.712 at the default); without the flux into the ODE, stable
  exactly for lam < pi^2, the claimed bound. The 2021b term format coupled the ODE to x(s = 0), zero
  by the BC, until 9bd434ab (Sachin, 12/04/2022) added D = 1 to match the batch format.
- Verdict: faithful to the data and interconnection of the source; its claim is contradicted at
  every lam ≥ 0, and the source does not fix the sign that would restore stability at the default.
  No edit (also maintainer, 11/19 04:00: PDE 8 left alone).

**X2**, `PIETOOLS_PDE_Example_Diagne.m` (Declan 2022-08-18 6a5d184e; Sachin 2022-11-22; Declan
2024-12).
- Source: Diagne, Bastin, Coron, "Lyapunov exponential stability of 1-D linear hyperbolic systems of
  balance laws", Automatica 48(1):109–114, 2012 (library [2]): d_t xi + Lambda d_x xi − M xi = 0,
  applied to the linearised Saint-Venant–Exner equations (abstract; full text PDF only).
- Comparison: the predecessor of `PIETOOLS_PDE_Ex_Diagne.m` (case 3, created by e21afbef,
  08/20/2022): same model and parameters. The header states the form [x_+(0); x_−(1)] in terms of
  [x_+(1); x_−(0)], but the code puts x1, of speed −3.76, in the x_+ slot, so the BCs prescribe each
  state at its outflow boundary (§3) and the PDE is not well posed forward in time. Case 3 prescribes
  the inflow boundaries since 598dfeb3 (Dan Braghini, 05/31/2026, permutation Ppm), gains k1 = −0.5,
  k2 = −4.5 (X2: 1, 2). No toolbox code calls X2; the 11/18 harness does (below).
- Verdict: unjustified (superseded duplicate, ill-posed). Moved with its GUI file to
  `cleaning/private/dead_code/PIETOOLS/PIETOOLS_examples/Examples_Library/Stability_Examples/` by
  `mv` (sha256 unchanged, index untouched); entry in `cleaning/README.md`.

References elsewhere (none edited). X2: `testing_codes/private/baselines/2026-11-18_studies/private/`
(`lb_meta.m` row X2; `st_run.m` line 261, `lb_catalog_add({'X1','X2','X3'})`, which now fails on X2;
the X2 rows of `results_NEWT*.csv`, `results_T1P*.csv`), and the records `PHASE_B.md`, `PHASE_C.md`,
`PHASE_C_TABLE.md`, `phaseC_plan*.csv`, `phaseC_plan.pl`, `phaseC_table.pl`,
`sol55/package_files.txt`, `settings_1116/psatz_faces_results/DECISIONS.md`. Cases 4.5, 8, 15, 16:
the library call sites, `lb_meta.m` rows 15 and 16 of the four harness folders, the command-line
scripts `stability_4.m`, `stability_8.m`, `stability_15.m`, `stability_16*.m`, and the PIESOS dead
code in `cleaning/private/dead_code/PIETOOLS/PIESOS/Examples/`. The manual cites cases 5 and 7 only.
Outside this check, the unused predecessors `PIETOOLS_PDE_Example_Saba.m`,
`PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m`, `PIETOOLS_PDE_Example_Lamare.m` and
`PIETOOLS_PDE_Example_Transport_Eq.m` (GUI files beside them; `stability_1.m` loads
`PIETOOLS_PDE_Example_Transport_Eq_GUI.mat`) stay for a maintainer decision.

PDFs not read: arXiv 1812.06794v1 (https://arxiv.org/pdf/1812.06794v1, the 12/2018 version of [8]);
Saba et al. 2020, Fig. 5 (https://doi.org/10.1109/TAC.2019.2934384); Diagne et al. 2012, the
Saint-Venant–Exner BCs and gains (https://perso.uclouvain.be/georges.bastin/paper95.pdf).

## 11. Removals, duplicates and the PDFs (11/19 18:30–19:00)

- **Case 16 removed from the stability examples** (maintainer: "remove it from the stability
  examples"). `PIETOOLS_PDE_Ex_Timoshenko_Beam_3.m` and its GUI file, with the command-line copies
  that carry the same BCs (`stability_16.m`, `stability_16_2.m`), moved to `orphaned_file/private/`
  (`mv`, checksums equal; the git index is unchanged). The library keeps `case 16` with an error
  that names the removal (checked: "Example 16 was removed from the stability examples on
  11/19/2026."). `stability_16_1.m`, `stability_16_3.m` and `stability_16_4.m`, with the beam BCs,
  stay; in `stability_16_4.m` the BC at s = 1 reads x4_sss − x3_s where its comment has gam·x3_s.
- **Predecessor duplicates moved to `orphaned_file/private/`** (maintainer: "move them to orphaned
  code if duplicates"), each compared with its library file, comments stripped:
  `PIETOOLS_PDE_Example_Saba.m` (case 4; differs by st, form, GUI paths and the typo sig2 = 2.048)
  with its five GUI files; `PIETOOLS_PDE_Example_Lamare.m` (case 2; differs by st and form) with its
  four GUI files; `PIETOOLS_PDE_Example_Transport_Eq.m` (case 1 mirrored: x_t = x_s with the BC at
  s = 1; its GUI file stays, `stability_1.m` loads it). `orphaned_file/private/README.md` lists them.
- **`PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m` stays**: not a duplicate. It lacks the derivative
  on the PDE term of the ODE (`term{2}.D = 1`, added in case 8 by 9bd434ab, 12/2022), so its ODE input
  is x(a), which the BC sets to 0: the ODE is decoupled, and the cascade is stable for lam < pi^2
  exactly, the documented bound. The flux input of case 8 (as its header documents) breaks that bound
  with either sign; with −Bxr (applied, §4) the threshold is 9.7946.
- **The PDFs** (maintainer: "you have my permission"), saved in the session scratchpad:
  - arXiv 1812.06794v1 (171 KB), §8.1.0.2: the mixed hyperbolic-diffusive Timoshenko form has
    A2 = [0;0;1;0] on u4_xx; it "does not appear to be stable in the given states" at c = 0, and the
    paper reports stability for any c > 0. The correction of case 15 is the system of the source; at
    c = 0 the 'Q' form now certifies it (§4), where the code of the paper found no certificate.
  - Diagne, Bastin, Coron, Automatica 48 (2012), §4.3, eq. (19) and the matrix K: each Riemann
    coordinate is prescribed at its inflow boundary (xi1, lam1 < 0, at x = L; xi2, xi3 at x = 0). X2
    prescribed them at the outflow boundaries; case 3 follows the paper. The paper gives no values
    for the speeds or the gains k1, k2: those of case 3 are its own.
  - Saba et al. 2020 (IEEE TAC): the IEEE PDF link returned nothing from this machine (no access);
    not fetched elsewhere. Case 4.5 rests on Peet, Automatica 2021, §10, Ex. 5, and on the computed
    threshold sig2 = 1.0484.

## 12. Notes removed; a recommended st in every example (11/19 19:15–20:00)

- **Settings notes removed from the example files and blocks** (maintainer: they "say nothing about the
  problem itself, just numerical tests and references to files the user doesn't have"). 57 header notes
  and 40 block notes (587 lines, `notes.pl`), the 11/17–11/19 settings paragraphs of six mixed-domain
  files (154 lines, `notes2.pl`), the "No second output st" paragraphs and the 11/18 st-block comments of
  21 files (`notes3.pl`). Comment lines only: the code lines of every file compare equal before and after.
  Where removed lines had been committed, the header log says so. Left: the dated test records of the
  posopvar session in its examples (other paragraphs of those headers) and older log entries that name
  `cleaning/private/` (DDE library 11/18, PDE library 11/18).
- **Base tiers of the st blocks before the additions:** example files stripped 26, light 24, heavy 6,
  veryheavy 4, extreme 3; DDE/NDS blocks stripped 31, light 3, heavy 1, veryheavy 3, extreme 2. The
  light-plus-deltas rewrite was withdrawn by the maintainer.
- **A recommended st in the 12 example files that had none** (maintainer: "Every example needs a
  recommended settings"; `addst.pl`, checked by `st_check12.m`: 12 of 12 return the intended st and form):
  - PDE 65 (Turing pair): light with `settings_2d_stock('light')`, form 'P': certified at 1,635,656
    variables, m 18,182, 387 s (STOCK2D; the posopvar session certified it so on 11/17 at 2,102,408).
  - PDE 55: stripped, 'P' (Sol job 2067511 pending).
  - PDE 36, 39, 41, 44, X3, 56, 69, 70: light, 'Q'; no measured setting certifies them. PDE 39 in 'P'
    with the stock storage (the posopvar session's untried candidate): unknown, ratio 0.498, 233,856.
  - The two Jagt-Peet files (`2D_Wave_Eq_Damped`, `2D_Plate_Eq`), which flag no test: light. The damped
    wave certifies in 'Q' at light (288,400 variables, 27 s, ratio 1.10); 'P' light is infeasible
    (754,180). The plate is not measured (one solve at d = 0 exceeded 17 minutes on 09/28).
  - The library collects st at the nine call sites (36, 39, 41, 44, 55, 56, 65, 69, 70).
- **Posopvar session reply (11/19 19:00):** PDE 65 was certified on 11/17 with the stock 2-D storage in
  'P'; PDE 39 on HEAD 33d30877 by the storage-based 'Q' removed from the executive on 11/17; nothing ever
  certified 56, 69, 70, 36, 41 or 44 in its library form.
- **PDE 65 on the library path** (LIBCHECK65, 19:33–19:41, `libcheck_fix3_cases.csv`): the file's own st
  builds the same program (1,635,656 variables, m 18,182); 'P' optimal, residual 6.6e-7, eta 2.8e-10,
  build 36 s, solve 427 s, 467 s in all.
- Two log lines of 11/19 (DDE library, NDS library) no longer carry gamma and variable counts.
- **Not removed (maintainer to decide):** 20 dated test-record paragraphs, 119 lines in 18 files, written
  by the posopvar session on 11/16–11/17 ("Certified (PIETOOLS_stability_sop, 'P' form, ...,
  MOSEK, 11/17/2026) ... (1,668,608 decision variables)"; references to `gap_results.md`,
  `S17_obstruction.md`, the settings proposal of 11/15). The problem facts in the paragraphs before them
  (thresholds, spectra, sources) would stay. Committed log entries that name `cleaning/private/` (DDE
  library 11/18, PDE library 11/18) are left as history.
- **Removed, 11/19 21:30** (maintainer: "Yes, strip them"; `notes4.pl`, explicit ranges checked against
  their first lines): the 20 paragraphs, 133 lines in 18 files with their blank separators, each file
  with a header log entry (every removed line was committed). Kept: the cross-references (the Riemann
  form of the damped wave with interior delay; library example 44 for the Xu system) and, reworded, the
  claim that no coercive PI operator P certifies the transport–edge-diffusion system in the form 'P' at
  any r2, marked "[argued, not checked by a second reader]", without its file reference. The st help
  line of the Laplace–Beltrami file no longer points at the header. `chk18.m`: 18 of 18 parse, resolve
  to the working copy and return the same PDE_t and st as the backup (scratchpad `examples_bak_notes4`).
  Left in example comments: the four committed log entries that name `cleaning/private/`.
