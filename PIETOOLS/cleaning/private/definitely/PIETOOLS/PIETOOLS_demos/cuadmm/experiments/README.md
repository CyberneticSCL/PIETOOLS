# experiments/

The one-off studies behind specific claims in `../BASELINE_REPORT.md` and in the campaign
notes. They are kept as the **record of how each result was measured**, not as a test suite
— for that, use `bl_check` and the lists in `../SUITES.md`.

Each script's header comment says what question it asks and what it found. By theme:

| theme | scripts |
|---|---|
| what SDP each executive builds | `census1d`, `census2d`, `shapemap1d`, `shapemap2d`, `probe_m`, `probe_psatz`, `suite_impact` |
| psatz generators and degree theory | `poinc1d`, `poinc2d`, `psatz1d`, `psatz2d`, `psatzgrid`, `psatzdeg`, `perdir`, `diag_deg`, `slot101`, `degmap`, `matchdeg`, `transfer2d`, `genplant`, `normgen`, `multprobe`, `odeprobe`, `odegen`, `localcoef`, `psz_blocks`, `psz_reach` |
| early validation of the instruments | `t0_path` … `t14_cvt`, `cuimport_test`, `twosided_export` |
| residual gates and ranking | `gatecmp`, `gatexfer`, `gaterank`, `gaterank2`, `gate_popcheck` |
| H∞ and bisection | `hinf_bisect`, `hinf_dump`, `cu_hinf_gam`, `hinf_ref`, `bis_dump`, `bl_fixtest`, `bl_h2probe`, `bl_h2ctrl` |
| nonlinear (PIESOS) | `fisher_dump`, `fisher_ladder`, `fisher_t1`, `fisher_time`, `run_fisher`, `fixednorm` |
| build walls and regressions | `dbg_build32`, `dbg2`, `hw_probe`, `hw_repro`, `regress_lpieq`, `verify_patch`, `verify_final` |
| cuADMM configuration | `dumpcfgs`, `dump1d`, `search2d`, `t13_mosek` |

## Running one

```matlab
cuadmm_path
run('experiments/poinc1d.m')
```

- They reach the package's `private/` helpers through `cuadmm_private('name', …)`.
- They write any generated data under `cuadmm_outdir()`, never into this folder.
- Scripts that **read** previously generated data — dumps under `baseline/dumps`,
  `cu_fish`, `cu_bis` and so on — expect it under `cuadmm_outdir()`. Those dumps were not
  kept (about 1.9 GB), so regenerate them first with the harness (`bl_run`, `bl_big`,
  `bl_bisdump`), or point `CUADMM_OUT` at a folder that still has them.
- Some scripts install a **shadow** (`cuadmm_shadow('lpisolve')` or `…('poslpivar')`).
  It stays active until `cuadmm_shadow('off')` or the next `cuadmm_path` — see the Shadows
  section of `../README.md`.

**None of these have been re-run since the move from the scratchpad.** The move rewrote their
paths only.

- **Most likely to work unchanged:** scripts that build their programs from scratch, such as
  `poinc1d`, `psz_reach` and `gatecmp`.
- **Won't run until their data is regenerated:** scripts that read old dumps.
- **Superseded:** a few of the psatz scripts predate the `poslpivar_2d` patch (commit
  `0cf1add1`). They were written against a shadow copy and may now exercise the committed
  version instead.
