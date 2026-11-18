# Mixed-domain diagnostics: S17 and S18 (2026-11-16)

Question (maintainer, 11/16/2026, after the psatz-faces run): on the two mixed-domain
systems that no setting has certified, is the obstruction the term set, the storage, the
degree, or the solver? S17: a 1-D transport state as Dirichlet data on the edge s2 = 0 of
a state diffusing in s2 only (order 0 in s1), margin r2* = pi^2. S18: library PDE 43 with
the sign of its header (a0 = 1.9, a1 = -1), the delay a 2-D transport state on
s2 in [0, tau], margin tau* = 1.034726. Both at half their margin: r2 = 0.5 pi^2,
tau = 0.517363. Light tier, dw pinned, executive path (PIETOOLS_stability_sop 'P'), the
scratch shadows of the psatz-faces run (harness/ of the parent folder) plus the files
listed below.

Runs (logs in the session scratchpad `psatz_faces/logs/m*.log`, records
`psatz_faces/results/exec_m*.mat`; scripts copied to `harness/` here):

- M1 `pf_run_m1.m`: term sets none, box (Putinar, per-direction offset 1), facesS
  (today's 'auto'), faces0, corners at the N-D default storage, dw 0 and 1, MOSEK.
- M2 `pf_run_m2.m`: SeDuMi on facesS and box at dw 1 (the plan's pairing rule for an
  UNKNOWN exit).
- M3 `pf_run_m3.m`: the storage with every basis operator of every space at the light
  degrees (PF_POS_SEP, route_sep shadow of exec_tools_sop): 'sep11' both directions
  separable (full integrals; the control for the extra mixed operators), 'sep01' s1
  split into lower and upper integrals, 'sep10' s2 split, 'sep00' both split; facesS,
  dw 0 and 1.
- CHK `pf_run_chk.m` / `pf_check_fam.m`: the certified derivative operators evaluated
  on the two-direction family of the obstruction argument (below).
- M6 `pf_run_m6.m`: S18 with every basis operator at raised storage degrees (mult 1, int 2 per
  variable, joint 4; PF_POS override), s2 split and both split.
- M7 `pf_run_m7.m` / `pf_box.m`: the stability LPI with ep I <= P <= I (ep 0.5, 0.1, 0.01) on
  S17, D4 and S3: the scale of P fixed, so a weakly infeasible program becomes Farkas-infeasible.
- M8 `pf_run_m8.m`: S18 attribution (separable storage at the raised degrees) and soundness
  (1.1 tau*, 0.9 tau*).
- M9 `pf_run_m9.m`: S18 with the two-operator default storage at the raised degrees.
- M5 `pf_run_m5.m` / `pf_dmax.m` (written, not run): the decay-margin program, maximize delta
  subject to P >= I and -(T'PA + A'PT) >= delta T'T.

Tables: `m1_m3_summary.md` (M1, M3, M6, M8, M9), `m2_chk_summary.md` (M2, CHK, CHK2, M7). Reading and the
obstruction argument: `S17_S18_findings.md`.

Second request (maintainer, 11/16 02:10): treat the failures as a parameterization defect
or the need for a weaker certificate, PIE2PDE being the default test; find fixed operators
whose positivity cannot be certified. Runs added (reading: `parameterization.md`):

- M10 `pf_run_m10.m`: the finite-energy form 'FP' on S17 (same exits as 'P').
- FIX `pf_run_fixed.m` / `pf_fixed.m`: fixed operators (T'T, T'E1T, T'E2T, A'A of S17, S3,
  S18; the certified S18 storage) under the slack reader: `fixed_operators.md`.
- Q `pf_run_q.m` with `qform/` (the executive with its 1-D-only restriction on the Q forms
  lifted, diff in `harness/qform/`): the storage-based Q form in 2-D.
- CONE `pf_run_cone.m` / `pf_cone.m`: fixed operators against the storage cones (X = P).
- QD `pf_run_qd.m` / `pf_qdirect.m`: the direct PIE2PDE form (no storage); QD2
  `pf_run_qd2.m`: its margin and large-margin soundness checks.

Logs: `psatz_faces/logs/{m10,next3,next4,qd,qd2}.log`; records `results/{fixed_mixed,
cone_mixed,qdirect_*,s18_Pstar}.mat`.

## 11/16, the term run and the tier baselines

Plan: terms_plan.md (+ Additions). Results: terms_results.md (T1: the lift of the PDE norm by term set; T2: the direct PIE2PDE test by term set), baselines_results.md (P2 (a)-(c), P3 (a)-(g): the candidate rows against the library, the mixed suite, the 2-D ladder, the 3-D heat, the per-direction cap, the heavy entry, the light candidates). Decisions D37-D53. Report: sopvar/settings_report_2026-11-16.pdf (its fragments settings_report_*.tex beside it; the construction fragment is re-extracted from sopvar_settings_proposal.tex by harness/fix_tiers_ref.pl). Harness: harness/ (pf_row_storage = the row reader; pf_rows = the candidate rows; pf_exec_row = one executive call under a row; the p2/p3 runners; the shadow diffs of poscopvar_direct (per-direction offset and cap) and exec_tools_sop (the row route)). Logs in the scratchpad only.
