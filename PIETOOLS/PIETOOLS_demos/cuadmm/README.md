# cuADMM evaluation for PIETOOLS

A harness for judging the GPU first-order SDP solver **cuADMM** against Mosek on the SDPs
that PIETOOLS executives generate, across every application class (stability, H∞,
estimator and controller synthesis, well-posedness, nonlinear local stability), plus a
settings ladder tuned to cuADMM's strengths. H2 was dropped from the testing regime on
2026-09-26; see `SUITES.md`.

- **`BASELINE_REPORT.md`** — what was measured and what it means. Start here.
- **`SUITES.md`** — the short test lists, and which to run after which kind of change.

**The baseline suite moved (11/08/2026).** The `bl_*` files except `bl_bisect`, `bl_face` and
`bl_fix` (which `heat2d_compare` and `jp_wave_plate_rates` call), `private/`, `results/`,
`SUITES.md`, `BASELINE_REPORT.md`, the `gpu/` shell scripts and `sol/` are in
`testing_codes/private/PIETOOLS_demos/cuadmm/` (maintainer: move the testing codes to
`testing_codes`; on 11/09/2026: "To be safe, move the files in testing folder a private
subfolder."). That folder is off the path: from this folder, add it after `cuadmm_path`, as in
`cuadmm_path; addpath(fullfile(pwd,'..','..','testing_codes','private','PIETOOLS_demos','cuadmm'))`,
then run the quick start below. `cuadmm_path` no longer asserts `bl_run` and `bl_check`. The files
this README names in those groups are there.

## Quick start

```matlab
cuadmm_path                 % PIETOOLS + solvers on the path; asserts nothing is shadowed
bl_check('smoke')           % 4 cases, ~1 min with startup: is the pipeline alive?
bl_check('structure')       % 16 cases, ~1.5 min: did a core change alter any assembled program?
```

Before committing a change to a core data structure, run at least `smoke` + `structure`.

Generated data — SDP dumps, cuADMM problem directories, logs, check output — goes under
**`cuadmm_outdir()`**, which is `<tempdir>/cuadmm` by default and **never** this folder. One
baseline run wrote ~1.9 GB of dumps.

```matlab
setenv('CUADMM_OUT','D:\cuadmm')   % put generated data somewhere durable
```

## Container path (11/07/2026)

Since 3df33dcc (10/29/2026), `convert` returns a `pie_struct` of `copvar` operators, which the
stock executives reject (`Top.dim`), and the opvar removal retires `opvar/`, `lpi_programming/`
and `executives/`. The harness now runs on the container path. It was run on 11/07 in an export
of 8fd693bd with those three folders deleted.

- The builders run the container executive of each stock name in `bl_cases`
  (`private/bl_sop`): `PIE2PDEstability` → `PIETOOLS_stability_sop(PIE,st,'Q')`,
  `Hinf_gain_2D` → `PIETOOLS_Hinf_gain_sop(PIE,st,'P')`, and so on.
- `bl_settings` sets `st.sop.max_raise = 0`, so a case is one SDP, as with the stock
  executives. With raises, the dumped SDP would depend on the solver's first verdict.
- `stab_mirror` builds the first pass of `PIETOOLS_stability_sop` form `'Q'` without solving.
  Its `At` and `b` equal the executive's exactly at heavy n = 1, 2, 4 and light n = 1.
- `cuadmm_path` asserts single resolution of `poscopvar`, `lpi_eq_sop`, `lpi_ineq_sop`,
  `lpivar_cdopvar` and `exec_tools_sop` in place of the stock LPI functions.
- `bl_run` and `bl_check` record a solve with no point (sossolve returns `x = []` on primal
  infeasibility) instead of erroring; `bl_check` banks it as `sentinel_nopoint`.
- `bl_run`'s `eta` leaves out rows of rounding residue, as sossolve's `info.eta` does. The
  container 2-D programs carry such rows: `stab2_rd` reads 5.2e-2 with them, 6.4e-10 without.

**The container executives pose different SDPs.** `results/2026-09-24/` and `BASELINE_REPORT.md`
hold stock-path values; `bl_expect.tsv` was re-recorded on 11/08/2026 (below). Measured 11/07
with MOSEK, stock light settings:

| case | stock m, blocks, result | container m, blocks, result |
|---|---|---|
| `stab_rd1` | 130, [10 15 6], certified | 126, [10 21 14], certified |
| `hinf_rd1` | 242, [10 31 13 1], gam 0.1826 | 195, [10 25 17 1], gam 0.18248164 |
| `hinfco_rd1` | 176, numerr 2 | 138, gam 0.1824806 |
| `hinfduco_rd1` | 168, gam 8251.9, numerr 2 | 126, primal infeasible, no point |
| `ctrl_rd1` | 181, rel_b 1.1e-3 | 141, primal infeasible, no point |
| `stab2_rd` | 3456, [8 424], rel_b 2.08 | 2968, 6 blocks (largest 324), certified |
| `hinf2_rd` | 14005, [8 869], gam 1.272 | 3945, 7 blocks (largest 334), gam 0.17172217 |
| `hinf2nc_rd` | – | 3470, numerr 2, gam 756.6 |
| `hinf2du_rd` | – | 3941, primal infeasible, no point |

`hinf2_rd`'s closed-form gain is 0.1711.

**Settings from `lpisettings_sop` (11/08/2026).** The stock `lpisettings` and
`settings_PIETOOLS_*` leave the MATLAB path, so `bl_settings`, `cuadmm_settings`, `bl_big`,
`bl_scale` and `heat2d_compare` take the tier of the same name from `lpisettings_sop`. Measured
on the 24 registry cases with MOSEK: the objective cases pose the same SDPs as in the table
above. The feasibility cases pose the tier's first pass, at its `dw_feas` and `eppos_stab`
1e-2: `stab_rd1` m 116, blocks [10 14 7]; `stab2_rd` m 2044, largest block 144; both
certified. `stab2_rd_psz90` (0.9 λ*) is primal infeasible at that pass; the tier's raise, which
`bl_settings` turns off, certifies it at the 11/07 program (m 2968).

`bl_expect.tsv` was re-recorded on these tiers (11/08/2026) for the 27 cases of the suites
`structure`, `twod`, `objective2d` and `scaling`; `stab2_rd_psz` and `stab2_rd_psz90` are
banked for the first time. A second `bl_check` pass gives 27 PASS, four of them known failures:
`hinfduco_rd1`, `ctrl_rd1` and `hinf2du_rd` return no point, `stab2_rd_psz90` the trivial point.
The H2 rows (dropped from the regime) and the `nl_fisher_*` rows (which do not build since the
PIESOS API change of 10/15/2026) keep their stock-path values. `hinf2nc_rd` banks as `pass` on
its residual (6.5e-6), although MOSEK returns numerr 2 and gam 756.6.

Not ported; each still calls the stock code:

- `experiments/` (45 of 73 scripts), moved on 11/07/2026 to
  `cleaning/private/definitely/PIETOOLS/PIETOOLS_demos/cuadmm/`, off the MATLAB path, with the
  helpers only they called: `cuadmm_private`, `cuadmm_shadow` with `shadows/` (instrumented
  copies of `lpisolve` and the stock `poslpivar`), and `private/lpi_resid`, `opnorm_pi`,
  `opgram`, `pgridval`, `gaussleg`, `legval` (Galerkin matrices of `opvar` operators), now in
  `private_files/` there;
- `cuadmm_linear_psatz` (no caller): it stops on the missing `eq_opts_psatz`, see below.

`heat2d_compare.m` runs on the container path since the low-rank package was ported (11/08/2026,
c3bfac42): its LPI is the build of `pielr_lpi('stability')`, the form `'P'` LPI of
`PIETOOLS_stability_sop` with the stock slack. Measured to assemble (nd 4: m 5721, blocks
[8 744]); not re-run end to end.

`stab2_rd_psz` and `stab2_rd_psz90` do not apply `pz`: the 2-D settings carry no `eq_opts_psatz`
since the master merge of 10/27/2026, and the container's default slack holds the four linear
face generators in 2-D anyway. `stab2_rd_psz` poses the same SDP as `stab2_rd` (measured).

`bl_bisect` and `cuimport` compute `eta` over every row, so on the container 2-D programs their
`eta` includes the residue rows (5.2e-2 at MOSEK's `stab2_rd` point, against the 1e-7
certification threshold).

`nl_fisher_*` stop in `PIESOS/@polyopvar/mtimes` ("Right-multiplication of polynomials is
currently not supported", from `innerprod.m:82`) with or without `opvar/`: a PIESOS change,
not this port.

## Layout

| path | contents |
|---|---|
| `cuadmm_path.m` | path setup; run first. Edit the SeDuMi/Mosek folders at its top for your machine |
| `cuadmm_settings.m` | `cuadmm_settings(tier,dim)`, the tiered settings ladder — same interface as `pielr_settings` |
| `cuadmm_linear_psatz.m` | turns on the four linear 2-D psatz generators on a settings struct |
| `bl_cases.m`, `bl_suites.m` | the case registry, and the question-organised lists over it |
| `bl_check.m` | runs a list against `bl_expect.tsv`; `bl_check(s,'bank')` records new cases |
| `bl_run.m`, `bl_scale.m`, `bl_big.m` | Mosek arm, size ladder, and the rungs Mosek cannot solve |
| `bl_bisect.m` | **certified bisection** on an objective-form dump, cuADMM or Mosek probes: F/I only from verified certificates, leans steer, threshold and cap scale with the bracket interval. Rule and calibration in its header; defaults in `cuadmm_settings(1).cuadmm.bisect` |
| `bl_fix.m`, `bl_bisdump.m` | pin γ to pose the feasibility question; write the bisection programs |
| `dump2cuadmm.m`, `cuimport.m`, `bl_verify.m` | SDP → cuADMM input; cuADMM certificate → scored against the ORIGINAL data |
| `bl_expect.tsv` | banked expectations — a package input, kept in the repo |
| `gpu/` | the cuADMM arm: `bl_gpu.sh`, `bl_bisgpu.sh`, `bl_gpu_post.sh` |
| `private/` | plant builders (`bl_b_*`, `fisher_prog`), the executive map `bl_sop`, `bl_settings`, `sdpshape`, `stab_mirror` |
| `experiments/` | moved 11/07/2026 to `cleaning/private/definitely/PIETOOLS/PIETOOLS_demos/cuadmm/experiments/`: the one-off studies behind specific claims in the report, with their README; their helpers moved with them (`cuadmm_private.m`, `cuadmm_shadow.m`, `shadows/`, six files of `private/`) |
| `results/2026-09-24/` | the measured tables and raw logs behind `BASELINE_REPORT.md` |
| `cuadmm_outdir.m` | output root |

## The GPU arm

cuADMM runs out of process under WSL on dumped problems:

```bash
bash gpu/bl_gpu.sh <dumps-dir> <out.tsv> 120 300
```

Two things that cost real time:

- It **needs `LD_LIBRARY_PATH=/usr/lib/wsl/lib`**; the scripts set it. A distro
  `libcuda.so.595.91.07` otherwise shadows the WSL driver, and every CUDA call reports "no
  CUDA-capable device".
- `cuadmm_exe` **exits 0 on total CUDA failure**. The scripts judge success by a
  `CUADMM_TIMING` line, never by the exit code.

The binary path defaults to this workstation's build; override it with `CUADMM_EXE`.

To build the binary elsewhere, run `SOSTOOLS400/internal/processing/cuadmm_install.sh`
(11/08/2026) under Ubuntu or WSL; on Windows, in PowerShell, in that folder:
`wsl bash cuadmm_install.sh`. It installs the missing packages, builds upstream cuADMM
(commit fcadc30) with `cuadmm_exe.patch` (11/07/2026) from the same folder into `~/cuADMM`,
where `cuadmm_sedumi` looks when `CUADMM_EXE` is unset, and solves a test SDP. The patch holds
the changes that `cuadmm_sedumi` and these scripts need: command-line arguments, the solution
dump, the `CUADMM_TIMING` line, and build fixes; its first lines give the manual build steps.
`gpu/cuadmm_pinf_inf.patch` applies on top of it.

The workstation build carries `gpu/cuadmm_pinf_inf.patch` (2026-09-27). It adds an optional
10th argument: `1` makes cuADMM's primal residual, both its stopping test and the printed pinf,
the worst row, ‖b−AX‖∞/(1+‖b‖∞). `bl_bisect(...,struct('pinf_norm','inf'))` uses it. The default
`0` reproduces the unpatched solver exactly (checked: `stab_rd1_hv` at 1e-6, 4023 iterations
before and after). Why: b is 88–99.8% zeros, so the 2-norm ratio is an absolute test that
tightens like 1/√m.

## Shadows: read this before running an experiment that uses one

`cuadmm_shadow`, `shadows/` and the experiments that used them moved to `cleaning/` on
11/07/2026. `cuadmm_path` still removes any shadow folder left on the path. The rest of this
section describes them as they were.

Some experiments replace a PIETOOLS function with an instrumented copy. `lpisolve` becomes a
capture stub that records the program and never solves; `poslpivar` gets a variant with
function-handle multipliers. `pietools_path_update` adds this whole repository with
`genpath`, so a real `lpisolve.m` anywhere in it would replace `lpisolve` **for every
PIETOOLS user**. The copies are therefore stored as `shadows/<name>.m.txt`.
`cuadmm_shadow(name)` writes one out as a real `.m` into a temp folder and puts that folder
first on the path, only for the experiment that asks for it.

A shadow **stays active until removed**. Every experiment that installs one now removes it
again at its end, but one that errors partway leaves it in place. While the `lpisolve` stub is
active, every call answers with a zero solution and a clean-looking "feasible" status. So the
stub now **warns on every call**, and an executive run by accident afterwards reports γ = 0
with a warning rather than silently.

`cuadmm_path` **removes** any active shadow, with a warning, every time it runs. That includes
the harness runners (`bl_run`, `bl_check`, …), which all call it. Never run a runner while you
want a shadow to stay active. Its single-resolution check can't see a shadow, because it
removes shadows first.

## What the move from the scratchpad changed (2026-09-25)

The harness was built and run in a flat, session-temporary scratchpad and moved here
afterwards. The move rewrote **paths only**:

| before | after | why |
|---|---|---|
| `pp;` (hardcoded paths, always `restoredefaultpath`) | `cuadmm_path;` | locates PIETOOLS itself; resetting the path is opt-in via `cuadmm_path('reset')` |
| `fullfile(HERE, …)` for generated data | `fullfile(cuadmm_outdir(), …)` | 1.9 GB of dumps must not land in a git tree |
| helper calls from `experiments/` | `cuadmm_private('name', …)` | `private/` is invisible from `experiments/` |
| `addpath(…'stub'…)`, `addpath(…'shadow'…)` | `cuadmm_shadow('lpisolve')`, `cuadmm_shadow('poslpivar')` | see Shadows |
| hardcoded `cuadmm_exe` | `${CUADMM_EXE:-…}` | per-machine binary |

**Not moved:**
- the ~1.9 GB of SDP dumps (they can be regenerated);
- the one-off `chain*.sh` orchestrators, which hardcoded scratchpad paths;
- a 12-line debug script.

**The harness has not been re-run in this location.** Every number in the report was measured
before the move. The first `bl_check('smoke')` here is also the check that the move broke
nothing.

## Status of the settings ladder

`cuadmm_settings` is **written but not validated**. Its header records what carries over from
the low-rank ladder. It also records that in 1-D the psatz is **necessary** with a searched
Lyapunov operator but **inert at `stripped`**, where its block is 1×1. That's measured by the
low-rank session's preset × flag cross, not re-run here. And it records an **inference from the code, not a measurement**: in 2-D both
ladders collapse to three distinct programs, because the 2-D executives discard every
top-level settings field. The reasoning, with line references, is in the header of
`cuadmm_settings.m`.

## Conventions

Files follow the sibling package `lowrank_2d_stability`: LF line endings (`.gitattributes`
pins them for `*.sh`, which break under WSL with CRLF), no licence block, and a
`% CC, date` marker on the function line of new entry points. The files moved from the
scratchpad are new to the repository, so they don't carry CLAUDE.md §1's per-line
modification markers. The move rewrote their paths only, and the table above records every
rewrite.
