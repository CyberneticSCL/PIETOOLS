# Example-specific settings (MMP, 11/18/2026)

Maintainer: "design settings files which are example specific, work better than the default tiers for
that application (e.g. lower cost for same accuracy), and are set inside the example file."

Written incrementally; every number below was measured on this workstation (MOSEK 11, MATLAB R2025b,
tree 33d30877 plus the uncommitted 11/17 tier work), with the library pass holding the other MATLAB slot.
Scripts and logs: `C:/Users/mpeet/AppData/Local/Temp/claude/363c499c6332a503/af90e0a1-bcff-4e8d-9a24-9aa4c52491cc/scratchpad/exst/`
(`exst_dde.m`, `exst_fix.m`, `exst_plan.m`, `exst_bisect.m`, `exst_exec.m`; logs in `logs/`).

## 1. Design: which knobs an example-specific setting carries

The library runs stability in the form `'Q'` of `PIETOOLS_stability_sop`. Form 'Q' is a generalized storage:
the storage is T'Q with Q free at `lpi.var.deg{nv}` (T'Q self-adjoint by an equality, P = (T'Q+Q'T)/2),
and the storage form `'P'` is its special case Q = P T. Q contains the P form of a tier when, per
direction, deg Q >= deg P + deg T, with deg P the kernel degree of the tier storage row (light: multiplier
and integral kernels of degree 2) and deg T the polynomial degree of the PIE operator T. The slacks of form
'Q' are sized from their own support. Read from the code
(`PIETOOLS_stability_sop.m`, the `'Q'` case and the `dwbase` block):

| knob | read in form 'Q'? | effect |
|---|---|---|
| `lpi.var.deg{nv}` | yes | degrees of Q: `[d1 d2 d3]` = multiplier, output and input degree d1, integral kernel degrees `[d3 d2]` (`lpivar_cdopvar`). The slack degrees follow from the support of T'Q, so this one knob sets the size of the whole program. |
| `sop.max_raise_feas(nv)` | yes | passes of the raise ladder after an INFEASIBLE/UNKNOWN pass; each raise rebuilds and resolves with the slack weight one higher. A failing parameter value costs `1 + max_raise_feas` solves. |
| `lpi.ineq.psatz{nv}`, `lpi.ineq.dD`, `lpi.ineq.dw_space` | yes | term set, lift raise, per-space weight raise of the slacks. |
| `sop.dw_feas`, `dw_P`, `dw_Q` | no | `dwbase = 0` in form 'Q'. |
| `lpi.pos` (storage rows), `settings_2d` | no | form 'Q' parametrizes the storage through Q, not through the storage row; the row matters only through the rule deg Q >= deg P + deg T. |

Measured consequence (before 01:24, when `lpi.var.deg` set Q): in 1-D, heavy and light with `lpi.var.deg{1} = [4 2 2]` build the same program at
every bisection step (equal decision variables and rows at every step, equal thresholds, Section 2), although
heavy differs from light also in storage rows, `dw_feas`, `dw_P` and `cap_delta`. An example-specific
setting therefore carries `lpi.var.deg{nv}` and `sop.max_raise_feas`, and, where a case needs it, the
slack options; it starts from `lpisettings_sop('light')`.

Gates. Harness gate: |feasratio-1| <= 0.3, no numerr, pinf, dinf (`lb_exec`). Executive gate: `info.status`
of the executive ('optimal' or 'inaccurate' certify) with the eta of the last pass. Both are given for
every certificate.

## 2. One-dimensional delay rows (DDE03, DDE04, DDE05, DDE07, DDE08)

Harness bisection (`ab_bisect` of the 11/17 A/B scripts with both gates per step; reltol 1e-3 of the
documented value, at most 14 solves; library bracket of the catalog). Variables and rows are those of the
last certified step (the threshold). "Solves" counts SDP solves including the raise passes; seconds are the
wall time of the whole bisection, one run each (the first variant of a session carries about 2-5 s of
warm-up). Every certificate listed passed both gates: harness |fr-1| <= 0.03 or better and executive
status 'optimal' with eta <= 4.4e-10.

Settings measured here through `lpi.var.deg{1}` and `sop.max_raise_feas(1)` (before the deg_Q separation of 01:24; Section 4). The settings finally set in the library blocks are in Section 7 (`sop.deg_Q{1}`, `sop.max_raise_Q(1)`).

| row | setting | threshold | variables | m | solves | seconds | fr / eta at threshold |
|---|---|---|---|---|---|---|---|
| DDE03 (doc 1.71785) | light | 1.6708775 | 982 | 155 | 23 | 9.3 | 0.999 / 2.2e-14 |
| | heavy | 1.7161724 | 1,712 | 179 | 15 | 4.4 | 0.989 / 4.4e-11 |
| | **[4 2 2], raise 0** | 1.7161724 | 1,712 | 179 | 11 | 2.7 | 0.989 / 4.4e-11 |
| | [3 2 2], raise 0 | 1.7144948 | 1,707 | 177 | 11 | 6.5 | 0.991 / 1.9e-12 |
| DDE04 (doc 6.17258) | light | 6.0580497 | 3,436 | 508 | 19 | 6.9 | 1.00 / 3.1e-11 |
| | heavy | 6.1665521 | 6,196 | 582 | 15 | 4.8 | 0.996 / 3.4e-12 |
| | **[3 2 2], raise 0** | 6.1665521 | 6,184 | 582 | 11 | 3.9 | 0.997 / 8.5e-12 |
| | [4 2 2], raise 0 | 6.1665521 | 6,196 | 582 | 11 | 3.4 | 0.996 / 3.4e-12 |
| DDE05 (doc 1.5707) | light | 1.5660983 | 859 | 128 | 17 | 4.7 | 0.995 / 3.6e-10 |
| | heavy | 1.5707 | 1,549 | 147 | 29 | 7.4 | 0.996 / 1.6e-10 |
| | **[4 2 2], raise 0** | 1.5707 | 1,549 | 147 | 11 | 2.4 | 0.996 / 1.6e-10 |
| | [3 2 2], raise 0 | 1.5691661 | 1,546 | 147 | 11 | 2.9 | 0.995 / 3.1e-11 |
| DDE07 (doc 1.3722) | light | 1.3478511 | 3,436 | 500 | 20 | 10.9 | 0.994 / 1.6e-13 |
| | heavy | 1.3721821 | 6,196 | 582 | 23 | 9.5 | 0.999 / 4.7e-11 |
| | **[4 2 2], raise 0** | 1.3721821 | 6,196 | 582 | 11 | 3.5 | 0.999 / 4.7e-11 |
| | [3 2 2], raise 0 | 1.369749 | 6,184 | 578 | 11 | 3.2 | 0.998 / 2.5e-12 |

Findings (measured):
- The threshold is set by the degree of Q alone: `[4 2 2]` reproduces heavy on all five rows with the
  same variables and rows at the threshold; `[3 2 2]` reaches heavy on DDE04 and falls 0.10-0.26% short
  on DDE03, DDE05, DDE07, DDE08.
- At `[3 2 2]` and `[4 2 2]` no raise pass ever turned a failed step into a certificate (92 failed steps
  over the five rows with raises 1 and 2, every one still failing after its raises), so the ladder only multiplies the cost of each failed step by
  3. Raise 0 cuts the solves of a bisection from 15-29 (heavy) to 11-12. Raises do certify when an
  integral-kernel degree of Q (d2 or d3) is 1 (`[3 2 1]`: 15 steps certified only on the second pass; `[3 1 1]`,
  `[4 1 1]` also lose threshold: DDE04 5.4010 and 5.4131 against light 6.0580).
- Against heavy the example setting builds the same program at each step; the gain is the solve count
  (11 against 15-29) and the bisection wall time (2.4-3.9 s against 4.4-9.5 s). Against light it trades
  1.6-1.8 times the variables at the threshold for a threshold 0.30% (DDE05) to 5.9%
  (DDE08) higher.
- NDS05 and NDS08: light already reaches heavy (NDS05 1.5784509 both) or exceeds it (NDS08 light
  0.027523096, heavy 0.027388574) in results_A.csv; no example setting needed, none designed.

### DDE08 (Example E; doc tauNmax in [3.9, 3.95], "requires veryheavy")

| setting | threshold | variables | m | solves | seconds | fr / eta |
|---|---|---|---|---|---|---|
| light | 1.3435913 | 1,258 | 218 | 30 | 7.1 | 0.994 / 3.4e-11 |
| heavy | 1.42229 | 2,068 | 252 | 26 | 6.2 | 0.992 / 6.4e-12 |
| [4 2 2], raise 0 | 1.42229 | 2,068 | 252 | 12 | 2.7 | 0.992 / 6.4e-12 |
| veryheavy | 3.8956787 | 50,214 | 854 | 33 | 33.2 | 0.992 / 1.6e-9 (3rd raise) |
| light, [5 2 2] | 1.42229 | 3,833 | 302 | | 7.9 | 0.995 / 1.5e-10 |
| light, [4 4 4] | 1.42229 | 7,164 | 426 | | 9.2 | 0.997 / 1.3e-10 |
| light, [6 2 2] | 3.179895 | 10,502 | 415 | | 11.2 | 0.98 / 3.0e-12 |
| light, [6 3 3] | 3.6146118 | 13,344 | 470 | | 9.6 | 0.994 / 9.1e-12 |
| light, [6 4 4] | 3.6820679 | 14,826 | 526 | 26 | 11.2 | 0.996 / 1.3e-12 |
| light, [6 6 6] | 3.8582031 | 36,738 | 780 | 27 | 20.7 | 0.999 / 3.3e-10 |

The multiplier degree d1 decides DDE08: d1 <= 5 stays at heavy's 1.42229, d1 = 6 moves to 3.18-3.86.
Near the threshold most steps end MOSEK UNKNOWN with numerr (harness 'inaccurate', executive 'unknown').

## 3. Two-variable examples, single solve at the library default parameters

One `PIETOOLS_stability_sop(PIE,st,'Q')` per (example, setting), the PIE of `lb_build_system` with the
library defaults; size gate 2e6 variables / m 45,000. Seconds are the wall time of the executive call and
vary by up to a factor 2 between identical programs (the library pass ran in the other slot; PDE 43 `heavy`
46.0 s and `[2 2 2], raise 0` 14.0 s build the same program), so variables and rows are the cost measure
and seconds are indicative. "passes" counts raise passes (1 = certified or failed without a raise).

Notation: `dD0` is `st.lpi.ineq.dD = 0` (slack lift degree at its minimum Dmin, light has Dmin + 1);
`deg2 = X` is `st.lpi.var.deg{2} = X`; `ps none` is `st.lpi.ineq.psatz{2} = 'none'` (the plain term with
factor 1 only; light has `faces1`, the plain term plus the 2nv face terms with factors (s_d - a_d) and
(b_d - s_d) one weight degree lower).

### 3.1 Examples that light certifies (PDE 33, 34, 35, 37, 40)

| PDE | setting | harness (fr) | executive (eta) | variables | m | passes | s |
|---|---|---|---|---|---|---|---|
| 33 Transport | light | C (0.974) | optimal (3.7e-10) | 52,660 | 2,552 | 1 | 3.2 |
| | dD0 | C (0.982) | optimal (1.9e-10) | 16,740 | 1,252 | 1 | 0.7 |
| | deg2 = 0 | C (1.14) | optimal (7.3e-10) | 4,169 | 549 | 1 | 0.3 |
| | deg2 = [1 0 0] | C (1.02) | optimal (8.6e-10) | 12,816 | 861 | 1 | 0.5 |
| | ps none | C (0.911) | optimal (1.9e-10) | 38,260 | 2,204 | 1 | 1.8 |
| | deg2 = 2 (heavy) | C (1.02) | optimal (5.9e-10) | 210,361 | 7,123 | 1 | 11.0 |
| 34 KISS | light | C (0.996) | optimal (2.0e-10) | 72,100 | 3,281 | 1 | 2.8 |
| | deg2 = [1 0 0] | C (0.996) | optimal (3.4e-10) | 72,016 | 3,063 | 1 | 2.4 |
| | dD0 | C (1.01) | optimal (4.4e-11) | 120,740 | 2,753 | 2 | 3.4 |
| | deg2 = 0 | C (0.997) | optimal (2.9e-10) | 234,009 | 4,623 | 2 | 10.7 |
| | ps none | F (-1) | infeasible | 291,700 | 6,653 | 2 | 22.7 |
| 35 Reaction-diffusion | light | C (1.01) | optimal (7.9e-10) | 72,100 | 3,281 | 1 | 2.2 |
| | dD0 | C (0.998) | optimal (1.1e-10) | 16,740 | 1,345 | 1 | 0.6 |
| | deg2 = 0 | C (1.05) | optimal (1.1e-9) | 37,449 | 2,325 | 1 | 1.4 |
| | ps none | U (2.01) | optimal (5.3e-10) | 57,700 | 3,281 | 1 | 2.9 |
| 37 Parabolic | light | C (1.00) | optimal (3.6e-10) | 72,100 | 3,265 | 1 | 1.9 |
| | dD0 | C (1.00) | optimal (2.7e-10) | 16,740 | 1,329 | 1 | 0.5 |
| | deg2 = 0 | C (0.991) | optimal (1.6e-10) | 72,009 | 3,049 | 2 | 2.1 |
| | ps none | F (-1) | infeasible | 291,700 | 6,637 | 2 | 25.4 |
| 40 Heat with ODE | light | C (1.00) | optimal (6.6e-10) | 73,463 | 3,356 | 1 | 3.1 |
| | dD0 | C (1.00) | optimal (2.2e-10) | 17,399 | 1,402 | 1 | 0.9 |
| | deg2 = 0 | C (1.00) | optimal (3.5e-10) | 73,366 | 3,111 | 2 | 3.2 |
| | ps none | C (1.00) | optimal (1.0e-10) | 58,383 | 3,356 | 1 | 2.1 |

(C certified, F infeasible, U unknown under the harness gate.) Combinations of dD0 with deg2 = 0 and
with ps none: Section 3.3 when measured.

### 3.2 Examples where light fails or the file already carries an st (PDE 38, 43, 57, 58, 66, 72)

`raise 0` is `st.sop.max_raise_feas = [2 0 0]`. "file st" is the st the example file returned on 11/17.

| PDE | setting | harness (fr) | executive (eta) | variables | m | passes | s |
|---|---|---|---|---|---|---|---|
| 38 Wave | light (results_A) | U (-0.121) | – | 245,200 | 13,167 | | 34.3 |
| | heavy (results_A) | C (0.793) | – | 744,164 | 30,999 | | 260.6 |
| | deg2 = 2, raise 0 | C (0.793) | optimal (1.3e-9) | 744,164 | 30,999 | 1 | 289.5 |
| | deg2 = [1 2 2], raise 0 | C (1.05) | optimal (1.1e-9) | 744,000 | 30,933 | 1 | 224.6 |
| | deg2 = 2, ps none, raise 0 | C (1.27) | optimal (1.7e-7) | 595,684 | 30,999 | 1 | 212.3 |
| | deg2 = [2 1 1], raise 0 | U (-0.252) | optimal (2.8e-9) | 245,284 | 13,203 | 1 | 34.3 |
| 43 Heat, interior delay | light | C (0.999) | optimal (7.8e-11) | 374,827 | 6,576 | 2 | 39.4 |
| | heavy = deg2 = 2 | C (0.999) | optimal (3.7e-10) | 215,188 | 8,686 | 1 | 14.0-46.0 |
| | deg2 = [2 1 1] | C (0.998) | optimal (1.1e-10) | 374,875 | 6,593 | 2 | 19.1 |
| | deg2 = 2, ps none | U (-0.222, numerr 2) | unknown | 816,213 | 15,583 | 2 | 132.4 |
| 57 Damped wave, delay | light | C (1.01) | optimal (2.3e-9) | 433,725 | 7,871 | 2 | 30.0 |
| | heavy = file st (deg2 = 2) | C (1.10) | optimal (6.0e-9) | 248,497 | 10,311 | 1 | 20.8 |
| | **deg2 = [2 1 1]** | C (0.961) | optimal (2.3e-9) | 70,465 | 4,381 | 1 | 4.7 |
| | deg2 = 2, ps none | C (0.962) | optimal (5.1e-10) | 189,853 | 10,311 | 1 | 15.0 |
| 58 Same, Riemann form | light | F (-1) | infeasible | 563,813 | 8,126 | 2 | 36.8 |
| | heavy = file st (deg2 = 2) | C (0.999) | optimal (5.5e-11) | 302,957 | 10,064 | 1 | 21.2 |
| | **deg2 = 2, ps none** | C (1.01) | optimal (7.4e-11) | 232,777 | 10,064 | 1 | 28.8 |
| | deg2 = [2 1 1] | F (-1) | infeasible | 563,890 | 8,203 | 2 | 36.1 |
| 66 Age-structured | light | C (1.01) | optimal (2.4e-10) | 56,548 | 2,775 | 1 | 3.9 |
| | heavy = file st | C (1.00) | optimal (1.4e-8) | 210,361 | 7,284 | 1 | 20.4 |
| | **dD0** | C (1.00) | optimal (9.4e-11) | 16,740 | 1,301 | 1 | 1.3 |
| | deg2 = 0 | F (-1) | infeasible | 72,009 | 3,060 | 2 | 2.9 |
| | ps none | F (-1) | infeasible | 233,380 | 5,685 | 2 | 19.7 |
| 72 Transport-edge diffusion | light | U (0.619, numerr 2) | unknown | 676,410 | 7,989 | 2 | 34.1 |
| | heavy = file st (deg2 = 2) | U (1.78) | optimal (1.1e-9) | 267,968 | 8,755 | 1 | 18.5 |
| | deg2 = [2 1 1] | U (0.595) | optimal (1.3e-9) | 676,458 | 8,018 | 2 | 36.5 |
| | deg2 = 2, ps none | F (-0.962, numerr 2) | unknown | 1,023,708 | 15,984 | 2 | 174.1 |

Findings (measured): in form 'Q' the heavy tier and light with deg2 = 2 build the same program (equal
variables and rows on PDE 33, 34, 35, 37, 38, 40, 43, 57, 58, 66, 72; heavy rows of 33-40 from results_A.csv). The file st of 57, 58, 66 and 72 therefore costs exactly the
heavy tier; of 43 (`dw_space = [1 0]`), measured in Section 3.3. The file st of 72 certifies under the
executive gate only (fr 1.78).

### 3.3 Combinations (PDE 33, 34, 35, 37, 40; measured 01:30, after the deg_Q separation)

| PDE | setting | harness (fr) | executive (eta) | variables | m | passes |
|---|---|---|---|---|---|---|
| 33 | **dD0, deg_Q{2} = 0** | C (0.972) | optimal (1.7e-9) | 409 | 111 | 1 |
| | dD0, deg_Q{2} = [1 0 0] | C (0.999) | optimal (2.1e-9) | 1,376 | 193 | 1 |
| | dD0, ps none | C (1.07) | optimal (4.5e-10) | 13,412 | 1,252 | 1 |
| 34 | dD0, deg_Q{2} = 0 | F (-1) | infeasible | 48,649 | 1,751 | 2 |
| | dD0, deg_Q{2} = [1 0 0] | C (0.996) | optimal (6.8e-12) | 120,656 | 2,535 | 2 |
| | dD0, ps none | F (-1) | infeasible | 67,492 | 2,753 | 2 |
| | dD0, deg_Q{2} = 2 | C (0.99) | optimal (5.7e-10) | 522,441 | 7,290 | 1 |
| 35 | **dD0, deg_Q{2} = 0** | C (0.984) | optimal (1.6e-11) | 8,009 | 815 | 1 |
| | dD0, deg_Q{2} = [1 0 0] | C (0.997) | optimal (7.0e-12) | 16,656 | 1,127 | 1 |
| | dD0, ps none | C (1.00) | optimal (1.1e-9) | 13,412 | 1,345 | 1 |
| 37 | dD0, deg_Q{2} = 0 | C (0.987) | optimal (1.1e-10) | 48,649 | 1,737 | 2 |
| | **dD0, deg_Q{2} = [1 0 0]** | C (1.01) | optimal (9.1e-11) | 16,656 | 1,117 | 1 |
| | dD0, ps none | F (-1) | infeasible | 67,492 | 2,737 | 2 |
| 40 | **dD0, deg_Q{2} = 0** | C (0.961) | optimal (6.2e-11) | 8,518 | 844 | 1 |
| | dD0, deg_Q{2} = [1 0 0] | C (1.02) | optimal (2.3e-10) | 17,315 | 1,181 | 1 |
| | dD0, ps none | C (0.993) | optimal (2.6e-10) | 13,743 | 1,402 | 1 |

Chosen (bold), set in the example files: 33 409 variables against light 52,660 (129 times fewer); 35
8,009 against 72,100 (9.0 times); 37 16,656 against 72,100 (4.3 times); 40 8,518 against 73,463 (8.6
times). All at the library default parameters; no threshold was bisected for these, so the margin to the
stability limit at the cheaper setting is not measured (PDE 34 and 37 document a limit: lam <= 2 pi^2 at
lam = 19, a <= pi^2/2 at a = 4).

**PDE 34 (KISS): nothing beats light.** Lowest measured: deg_Q{2} = [1 0 0] at 72,016 variables against
light 72,100 (0.1% fewer, m 3,063 against 3,281); every cheaper slack (dD0, ps none, deg_Q{2} = 0) fails
the first pass. No st set.

## 4. The deg_Q separation of 11/18/2026 (01:24) and what it changes here

Since 01:24 the direct form reads its Q degree from `st.sop.deg_Q` and its raises from
`st.sop.max_raise_Q` (`raise_Q` 'slack', 'deg' or 'both'); `lpi.var.deg` is the synthesis degree only.
Every measurement above that set `lpi.var.deg` or `sop.max_raise_feas` was taken either before 01:24
(then those fields drove form 'Q') or, from 01:27 on, with the scratch helper setting `sop.deg_Q` and
`sop.max_raise_Q` alongside, so each variant means what its name says. One run straddled the change
(PDE 39, 01:23-01:28) and was discarded. The example files set `st.sop.deg_Q` and `st.sop.max_raise_Q`.
Since 01:29 the light tier carries 1-D deg_Q = [4 2 2] (maintainer); the example settings of the DDE
rows then differ from light only in `max_raise_Q(1) = 0` (and DDE04 in deg_Q{1} = [3 2 2]).

## 5. Form 'P' against form 'Q' on the 1-D delay rows (requested by the maintainer, 11/18)

Same harness bisection, `PIETOOLS_stability_sop(PIE,st,'P')` (the storage of the tier storage row),
light and heavy tiers as they stand at 01:33 (form 'P' does not read deg_Q).

| row | form 'P', light: threshold / variables / m | form 'P', heavy | form 'Q', smallest deg_Q{1} reaching the same threshold (variables / m) |
|---|---|---|---|
| DDE03 | **1.71785** / 1,220 / 70 (2nd raise) | 1.71785 / 1,245 / 70 | not reached at [4 2 2] (1.7161724, 1,712 / 179); Section 5.1 |
| DDE04 | 6.1665521 / 2,604 / 186 | 6.1665521 / 4,520 / 238 | [3 2 2]: 6,184 / 582 |
| DDE05 | 1.5707 / 651 / 47 | 1.5707 / 1,130 / 60 | [4 2 2]: 1,549 / 147 ([3 2 2]: 1.5691661) |
| DDE07 | 1.3721821 / 1,356 / 140 | 1.3721821 / 4,520 / 239 | [4 2 2]: 6,196 / 582 ([3 2 2]: 1.369749) |
| DDE08 | 1.42229 / 936 / 77 | **3.5471558** / 2,981 / 125 | light P: [4 2 2], 2,068 / 252; heavy P: [6 3 3], 3.6146118 at 13,344 / 470 ([6 2 2]: 3.179895) |

Every certificate in this table: harness fr in [0.998, 1.00], executive 'optimal', eta <= 1.7e-10.

Measured consequence: on all five rows the light tier in form 'P' reaches the threshold of the best form-'Q'
setting at 1.4-4.6 times fewer variables and 2.6-4.2 times fewer rows (variables of the form-Q setting over those of light P), and on DDE03 a higher threshold
(1.71785, the documented value). The form-'Q' family contains form 'P' only from deg_Q >= deg P + deg T on;
on these rows that inclusion costs more than the P form itself. An example setting cannot select the form
(the executive call names it: `PIETOOLS_DDE.m`, `PIETOOLS_auto_execute_sop`), so this is reported, not set.

## 6. Mixed-domain reruns after the separation (01:40-01:45) and form 'P' at light in two variables

| PDE | setting | harness (fr) | executive (eta) | variables | m | passes |
|---|---|---|---|---|---|---|
| 43 | 11/17 file st (dw_space [1 0]) | C (1.16) | optimal (4.8e-10) | 397,870 | 6,734 | 2 |
| | dD0 | F (-1) | infeasible | 115,647 | 3,282 | 2 |
| | dD0, dw_space [1 0] | F (-1) | infeasible | 124,954 | 3,380 | 2 |
| | **dD0, deg_Q{2} = 2** | C (0.998) | optimal (6.3e-10) | 196,513 | 6,376 | 1 |
| 57 | dD0 | C (1.03) | optimal (1.2e-8) | 140,457 | 4,109 | 2 |
| | **dD0, deg_Q{2} = [2 1 1]** | C (1.01) | optimal (1.7e-9) | 23,121 | 2,344 | 1 |
| | dD0, deg_Q{2} = 2 | C (1.00) | optimal (2.0e-8) | 240,830 | 7,821 | 1 |
| | deg_Q{2} = [2 1 1], ps none | C (0.793) | optimal (2.5e-9) | 241,722 | 7,953 | 2 |
| 58 | dD0 | F (-1) | infeasible | 201,325 | 4,246 | 2 |
| | dD0, deg_Q{2} = [2 1 1] | F (-1) | infeasible | 201,402 | 4,323 | 2 |
| | **dD0, deg_Q{2} = 2** | C (1.00) | optimal (2.4e-11) | 225,705 | 6,608 | 1 |
| | deg_Q{2} = [2 1 1], ps none | F (-1) | infeasible | 316,530 | 8,203 | 2 |
| 66 | 11/17 file st in form 'P' | C (1.00) | optimal (2.0e-10) | 189,216 | 3,022 | 2 |
| | heavy, dw_feas(2) 0, dD0, form 'P' | F (-1) | infeasible | 149,776 | 2,153 | 3 |
| | same, form 'Q' (= light, dD0) | C (1.00) | optimal (9.4e-11) | 16,740 | 1,301 | 1 |

Light tier in form 'P' (one solve, library defaults):

| PDE | deg T (s1, s2) | harness (fr) | executive | variables | m | passes | example st in 'Q' |
|---|---|---|---|---|---|---|---|
| 33 | 0, 0 | C (1.02) | optimal (9.9e-11) | 26,545 | 1,405 | 1 | 409 |
| 34 | 2, 2 | C (0.975) | optimal (1.5e-10) | 188,545 | 2,987 | 2 | none (light 'Q' 72,100) |
| 35 | 2, 2 | C (0.995) | optimal (1.5e-9) | 26,545 | 1,443 | 1 | 8,009 |
| 37 | 1, 1 | C (0.999) | optimal (1.5e-9) | 188,545 | 2,984 | 2 | 16,656 |
| 40 | 1, 1 | U (0.654) | optimal (2.1e-10) | 27,177 | 1,442 | 1 | 8,518 |
| 43 | 2, 0 | F (-1) | infeasible | 428,185 | 4,626 | 2 | 196,513 |
| 57 | 2, 0 | F (-1) | infeasible | 496,381 | 5,422 | 2 | 23,121 |
| 58 | 0, 0 | F (-1) | infeasible | 638,505 | 5,663 | 2 | 225,705 |
| 66 | 1, 2 | – | – | – | – | – | 16,740 |

deg T: the largest exponent of s_k plus that of s'_k over the monomials with a nonzero coefficient in
every block and multi-index of `PIE.T` (scratch `exst_degT.m`, not checked against a second reader).
In two variables form 'Q' wins against form 'P' at light on every case measured: the same verdict at fewer
variables (33, 34 at light 'Q', 35, 37, 40) or a certificate where 'P' at light is infeasible (43, 57, 58).
In one variable form 'P' wins (Section 5). The smallest deg_Q{2} that reaches light-'P' verdicts: 0 for
33, 35, 40 (with dD0), [1 0 0] for 37 (with dD0), the light value [1 1 1] for 34. These are below the
inclusion bound deg P + deg T (4 for 35 and 34, 3 for 37 and 40, 2 for 33), so on these cases the bound is
sufficient for inclusion but far from necessary for the verdict.

### 5.1 (Section 5 continued) Form 'Q' degrees between [4 2 2] and the light-'P' threshold (01:45, deg_Q, raise 0)

| row | [5 2 2] | [4 3 3] | [6 2 2] |
|---|---|---|---|
| DDE03 | **1.71785**, 3,337 / 221 | 1.71785, 3,370 / 253 | 1.71785, 4,174 / 244 |
| DDE04 | 6.1665521, 12,416 / 732 | 6.1665521, 12,540 / 866 | 6.1665521, 15,628 / 812 |
| DDE05 | 1.5707, 3,104 / 185 | 1.5707, 3,135 / 218 | 1.5707, 3,907 / 205 |
| DDE07 | 1.3721821, 12,416 / 732 | 1.3721821, 12,540 / 858 | 1.3721821, 15,628 / 812 |
| DDE08 | 1.42229, 3,833 / 302 | 1.42229, 3,870 / 332 | 1.42229, 4,738 / 331 |

(threshold, variables / m at the threshold; every certificate harness fr in [0.979, 1.00], executive
'optimal', eta <= 5.6e-10.) DDE03 needs deg_Q{1} = [5 2 2] for the documented 1.71785, which veryheavy
reaches at 63,810 variables (results_A.csv) and light in form 'P' at 1,220. DDE03 example st: [5 2 2]. On
the other rows degrees above [4 2 2] add variables and no threshold.

Rerun of the chosen settings through `sop.deg_Q` and `sop.max_raise_Q` (part 'v', 01:45): every
threshold, variable count and row count equal to the pre-separation runs of Section 2.

deg T on the DDE PIEs: the scratch routine returns 0 in s1 for all five rows; the maintainer states degree 1.
The routine reads the exponents from `ZL`/`ZR` of each block and may miss an exponent carried elsewhere;
this disagreement is unresolved, and the deg T values of Section 6 carry the same caveat.

## 7. Summary: the settings set in the example files (form 'Q', library defaults or harness bisection)

| example | st (from `lpisettings_sop('light')`) | st: result, variables / m | best tier: result, variables / m |
|---|---|---|---|
| DDE03 (library block) | deg_Q{1} [5 2 2], max_raise_Q(1) 0 | 1.71785, 3,337 / 221, 11 solves | veryheavy 1.71785, 63,810 (heavy 1.7161724, 1,712 / 179) |
| DDE04 | deg_Q{1} [3 2 2], max_raise_Q(1) 0 | 6.1665521, 6,184 / 582, 11 solves | heavy 6.1665521, 6,196 / 582, 15 solves |
| DDE05 | deg_Q{1} [4 2 2], max_raise_Q(1) 0 | 1.5707, 1,549 / 147, 11 solves | heavy 1.5707, 1,549 / 147, 29 solves |
| DDE07 | deg_Q{1} [4 2 2], max_raise_Q(1) 0 | 1.3721821, 6,196 / 582, 11 solves | heavy 1.3721821, 6,196 / 582, 23 solves |
| DDE08 | deg_Q{1} [8 6 6] | 3.9144165, 30,356 / 756, 22 solves, 17.1 s | veryheavy 3.8956787, 50,214 / 854, 33 solves, 33.2 s |
| PDE 33 | dD 0, deg_Q{2} 0 | certified, 409 / 111 | light certified, 52,660 / 2,552 |
| PDE 35 | dD 0, deg_Q{2} 0 | certified, 8,009 / 815 | light certified, 72,100 / 3,281 |
| PDE 37 | dD 0, deg_Q{2} [1 0 0] | certified, 16,656 / 1,117 | light certified, 72,100 / 3,265 |
| PDE 40 | dD 0, deg_Q{2} 0 | certified, 8,518 / 844 | light certified, 73,463 / 3,356 |
| PDE 38 | deg_Q{2} 2, psatz{2} 'none', max_raise_Q(2) 0 | certified (fr 1.27, eta 1.7e-7), 595,684 / 30,999 | heavy certified, 744,164 / 30,999; light UNKNOWN |
| PDE 43 | dD 0, deg_Q{2} 2 | certified, 196,513 / 6,376 | heavy certified, 215,188 / 8,686 |
| PDE 57 | dD 0, deg_Q{2} [2 1 1] | certified, 23,121 / 2,344 | heavy certified, 248,497 / 10,311 |
| PDE 58 | dD 0, deg_Q{2} 2 | certified, 225,705 / 6,608 | heavy certified, 302,957 / 10,064; light infeasible |
| PDE 66 | dD 0 | certified, 16,740 / 1,301 | light certified, 56,548 / 2,775 |
| PDE 72 | dD 0, deg_Q{2} 2 | executive 'optimal' only (fr 0.675), 241,648 / 6,257 | heavy executive 'optimal' only (fr 1.78), 267,968; light UNKNOWN |

The 11/17 st of 57, 58 and 72 set `lpi.var.deg{2}`, which form 'Q' no longer reads since the separation
(01:24); without the 11/18 change those examples would have run at light's Q degree.

Nothing beat the tier: PDE 34 (KISS; Section 3.3); the second gate on PDE 72 (no setting measured gives
|fr-1| <= 0.3); NDS05, NDS08 (light already at or above heavy). Not reached in the time: PDE 36, 39 (Section
8), the size-gated mixed cases 50, 61, 65, 69, 70, and PDE 45-48.

## 8. PDE 39 (2-D Navier-Stokes, Antonelli): no clean certificate found

deg T = 2 in s1 and s2 (scratch routine), so the inclusion bound for light's storage is deg_Q{2} >= 4 in
each direction. Measured, form 'Q', library defaults, raise 0 in two variables (02:00-02:10):

| setting | harness (fr) | executive (eta) | variables | m | s |
|---|---|---|---|---|---|
| light (results_A.csv) | U (7.14) | optimal | 648,900 | 30,721 | 260.5 |
| deg_Q{2} = [2 1 1] | U (2.65) | optimal (1.4e-8) | 649,089 | 30,725 | 297.8 |
| psatz{2} = 'none' | F (-0.999, numerr 2) | unknown (4.9e-9) | 519,300 | 30,721 | 335.2 |
| deg_Q{2} = [1 2 2] | refused by the size gate | – | 3,736,080 | 83,704 | – |

Raising the multiplier degree of Q to 2 moves OPT_STATUS from 7.14 to 2.65 at 189 more variables; raising
the integral-kernel degrees to 2 costs 3.74M variables, over the gate (degree 2 throughout: 3.7M, measured 11/17 per the brief). The storage form of HEAD (933,696 variables, certified) cannot be selected through st in the
direct form. No st set for PDE 39.

## 9. Run record

A queued chain of mine that I believed dead reran mixB (01:31), p38 (02:05) and mixC (02:17); I stopped
each. Their reruns overwrote the scratch logs `exst_fix_mixB.log`, `exst_fix_p38.log`, `exst_fix_mixC.log`;
the numbers of Sections 3.2 for PDE 38, 57, 58, 72 were copied from the first runs before that. One p39
run (01:23-01:28) straddled the deg_Q separation and was discarded (`exst_fix_p39_killed.log`). Runs between
01:25 and 01:29 errored on a parse error of `exec_sop_settings.m` (fixed by the caller at 01:29) and were
rerun.

## 10. Verification through the example files (02:00-02:32, tree dd465ee5 + working tree)

`[PDE_t,st] = PIETOOLS_PDE_Ex_<name>(0,{})`, `initialize_PIETOOLS_PDE`, `convert(PDE,'pie')`,
`PIETOOLS_stability_sop(PIE,st,'Q')` (scratch `exst_verify.m`): every file reproduces the table of Section 7
(variables, rows and both gates identical): 2D_Transport_Eq 409, 2D_Reaction_Diffusion_Eq 8,009,
2D_Parabolic_Eq 16,656, 2D_Heat_Eq_with_ODE 8,518, Age_Structured 16,740, Damped_Wave_Interior_Delay
23,121, its Riemann form 225,705, Heat_Eq_w_Interior_Delay 196,513, Transport_Edge_Diffusion 241,648
(executive 'optimal', harness fr 0.675), 2D_Wave_Eq 595,684 (fr 1.27, eta 1.7e-7). The DDE settings were
rerun through `sop.deg_Q`/`sop.max_raise_Q` in the harness bisection (part 'v', 'h'); the uncommented
library blocks themselves were not executed through `PIETOOLS_DDE.m`.

A last run (PDE 38 with dD0 and deg_Q{2} = 2; size-gated PDE 61 and 65 with dD0) was queued at 02:32. It never got a MATLAB slot, so I cancelled it at 02:44. Those cases were not measured.

## Maintainer decision, 11/18 02:55 (D100)
The form is not a settings knob: 'P' and 'Q' certify different forms of stability, and the 'P'
form, the stronger one, is to be used whenever it certifies. For the 1-D delay rows measured
here, where light in form 'P' reaches the thresholds at 1.4-4.6x fewer variables than any form-'Q'
setting, the example files should therefore name the 'P' form as their executive option rather
than carry a form-'Q' st; the 2-D blocks keep form 'Q' where 'P' does not certify (43, 57, 58).
