# Mixed-domain executive runs, 11/12/2026

The executive runs behind Layer 3 of
`testing_codes/private/converters/tests/mixed_domain/MIXED_DOMAIN_REPORT.md`:
the run of record of `test_mixed_domain_pde('all')` and every further run –
other parameters, tiers and corner rows, the subsystems of S17 and S18, the
control C1 and the runs near each margin.

**Trees.** `t`: `git archive` of `ndopvar` at 6c6f60dd, without
`testing_codes/`, `cleaning/`, `orphaned_file/`, `PIETOOLS_demos/` and
`PIETOOLS_documentation/`; the test folder of the live checkout is added to
the path. `t2`: `t` plus three uncommitted files of the session *New Posopvar
and Duality*, copied 03:40 (set A: `poscopvar_direct.m` cce8a04c0a,
`get_lift_degs.m` 5fc41d4f87, `lpi_ineq_sop.m` 5fd6f1dff8). `t3`: `t2` with
`get_lift_degs.m` 2cf69eda27 (set B, copied 04:34). On `t2` and `t3` the
runner sets `st.lpi.ineq.psatz = 'corners'`. `mx_frozen_path.m` sets the path
and asserts that every probed name resolves once, inside the tree.

**Machine.** Intel i9-14900KF (24 cores, 32 threads), 63.7 GB, Windows 11,
MATLAB R2025b, MOSEK 11.0. MOSEK ran at 8 threads: every log reports "TBB
Warning: The number of workers is currently limited to 7".

## Files

| File | Content |
|---|---|
| `results.csv` | One row per executive call of the runs N3–N18 (columns below): 73 rows, 68 returned; 2 stopped at a time cap (N9, N15b), 1 MATLAB crash and 2 runner errors (N17), each rerun or recorded as incomplete. |
| `test_runs.csv` | The check lines of the test: the run of record of `test_mixed_domain_pde('all')` (N7, 125 checks) and the C1 run (N13a, 20 checks); columns run, part, name, status, detail. |
| `private/` | The runners (`mx_*.m`) and queue scripts (`*.sh`) as they ran; the headers of the `.m` files were rewritten to the repository format after the runs, the code is unchanged. `mx_rows2csv.awk` and `mx_checks2csv.awk` wrote the two CSV files from the logs. The paths point at the session scratch folder: to rerun, export 6c6f60dd and change the root in `mx_frozen_path.m` and the scripts. |
| `logs/` | The MATLAB logs (not in git). |

## Runs

| Run | Script | Tree, tier | Content | Time (11/12) |
|---|---|---|---|---|
| N1, N2 | `night_now.sh` | t, light | fast layers and PIESIM; the exec layer, which crashed after 11 cases (R2025b access violation in `evalc`, since removed from the test) | 03:14–03:35; superseded by N7 |
| N3 | `night_now.sh` | t, light | S3, S8, S17 at 0.5 of the margin | to 03:40 |
| N4 | `night_now.sh` | t2, light, corner rows A | S3, S8, S17 at 0.5 and 0.9 | 03:40–04:30 |
| N5 | `night_now.sh` | t, heavy | S3 at 0.9 r*; S18 at tau = 0.5 | 04:30–05:22 |
| N6 | `night_now.sh` | t, veryheavy | S3 at 0.9 r* | 05:22–05:44 |
| N7 | `final_run.sh` | t, light | run of record: `test_mixed_domain_pde('all')`, 125 checks | 05:44–06:40 |
| N8 | `reach.sh` | t, light | S8 bisection over [0.5, 0.9] r3*, 4 steps; S17 at 0.25, 0.1 r2*; S18 at tau = 0.1, 0.25; S32 (S18 by the delay syntax), both signs | 06:40–06:55 |
| N9 | `reach.sh` | t, light | S13 (1-D -> 2-D -> 3-D chain): killed at the 90-min cap, no result | 06:55–08:25 |
| N10 | `corners2.sh` | t3, light, corner rows B | S3, S8 at 0.9 | 08:25–08:30 |
| N11a | `n11v3.sh` | t, light | C1: H-infinity 'Q', 'P'; H2 'o', 'oco' | 08:31–08:33 |
| N12 | `n11v3.sh` | t, light | D1–D5 (parts of S17, the S18 mode as a scalar delay equation); S18 with a0 = 0, a1 = 0.5 | 08:33–08:36 |
| N11b | `n11v3.sh` | t, heavy | S17, S8 at 0.9 | 08:36–10:06 |
| N11c | `n11v3.sh` | t, heavy | S24, H2 'oco' | 10:06–10:11 |
| N13a | `n13.sh` | t, light | C1 through the test (`convert`, `ops`, `exec`) and the whole `exact` layer | 10:11–10:14 |
| N13b | `n13.sh` | t, light | D6 (S18 x1 without the delay term); S18 with a1 = 0; S18 with a1 = 1 at a0 = 0.9, 1.2, 1.5 | 10:14–10:18 |
| N13c | `n13.sh` | t, heavy | C1: H-infinity 'Q'; H2 'o', 'oco' | 10:18–10:29 |
| N14 | `n14.sh` | t, light | C1: H2 'c' (the default form) and 'cco' | 10:29–10:33 |
| N15a | `n15.sh` | t, light | the near-margin set; MATLAB access violation after the first case | 10:33–10:34 |
| N15b | `n15.sh` | t, heavy | S3 at 1.02 r*: infeasible at dw 0 and 1, killed at the 40-min cap in dw 2; 1.1 r* not reached | 10:34–11:14 |
| N16 | `n16.sh` | – | skipped: N11b holds S8 at heavy | 11:14 |
| N17 | `n17.sh` | t, light | the near-margin set of N15a again, one MATLAB process per case: S26 absent from `mx_sys.m` (2 ERROR rows), S7 at 1.02 crashed (access violation) | 11:15–11:35 |
| N18 | `n18.sh` | t, light | S26 at 0.98 and 1.02 of its margin and S7 at 1.02, through the test catalogue (`mx_exec4.m`) | 11:36–11:47 |

Two starts overlapped and were killed: a second `night.sh` (identical to
`night_now.sh` but waiting for 04:30) began at 04:30, overwrote
`night_status.txt` and `N1.log`, and was stopped; an old queue script that
survived TaskStop started N11a a second time at 08:30, and both copies were
killed before N11a was rerun alone (08:31).

## Columns of `results.csv`

`run, tree, tier, system, params, executive, form, status, certified, dw,
ndv, m, seconds, gam, expected, gap, verdict, Ks, message`

- `params`: the system parameters (JSON), as passed to the system builder.
- `executive`: `stab` (`PIETOOLS_stability_sop`), `hinf`
  (`PIETOOLS_Hinf_gain_sop`), `h2` (`PIETOOLS_H2_norm_sop`); `form` the form
  passed.
- `status`, `certified`, `dw`: of the last solve; the executive raises dw
  while the solver returns UNKNOWN or INFEASIBLE, up to the limit of the tier.
- `ndv`, `m`, `Ks`: decision variables, equality rows and the first six PSD
  block sizes of the last solve; `seconds`: the executive call.
- `expected`: for `stab`, 1 if the PIE is stable and 0 if not; otherwise
  the exact gain or H2 norm. `gap` = gam/expected - 1.
- `verdict`: `pass` (stable and certified, unstable and not certified, or a
  certified bound at or above the exact value), `conservative` (stable, not
  certified), `not certified` (gain or norm), `UNSOUND` (an unstable PIE
  certified or a certified bound below the exact value), `ERROR`.
