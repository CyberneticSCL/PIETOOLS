# cuADMM evaluation for PIETOOLS

A harness for judging the GPU first-order SDP solver **cuADMM** against Mosek on the SDPs
that PIETOOLS executives generate, across every application class (stability, H∞, H2,
estimator and controller synthesis, well-posedness, nonlinear local stability), plus a
settings ladder tuned to cuADMM's strengths.

- **`BASELINE_REPORT.md`** — what was measured and what it means. Start here.
- **`SUITES.md`** — the short test lists, and which to run after which kind of change.

## Quick start

```matlab
cuadmm_path                 % PIETOOLS + solvers on the path; asserts nothing is shadowed
bl_check('smoke')           % 4 cases, ~1 min with startup: is the pipeline alive?
bl_check('structure')       % 22 cases, ~1.5 min: did a core change alter any assembled program?
```

Before committing a change to a core data structure, run at least `smoke` + `structure`.

Generated data — SDP dumps, cuADMM problem directories, logs, check output — goes under
**`cuadmm_outdir()`**, which is `<tempdir>/cuadmm` by default and **never** this folder. One
baseline run wrote ~1.9 GB of dumps.

```matlab
setenv('CUADMM_OUT','D:\cuadmm')   % put generated data somewhere durable
```

## Layout

| path | contents |
|---|---|
| `cuadmm_path.m` | path setup; run first. Edit the SeDuMi/Mosek folders at its top for your machine |
| `cuadmm_settings.m` | `cuadmm_settings(tier,dim)`, the tiered settings ladder — same interface as `pielr_settings` |
| `cuadmm_linear_psatz.m` | turns on the four linear 2-D psatz generators on a settings struct |
| `bl_cases.m`, `bl_suites.m` | the case registry, and the question-organised lists over it |
| `bl_check.m` | runs a list against `bl_expect.tsv`; `bl_check(s,'bank')` records new cases |
| `bl_run.m`, `bl_scale.m`, `bl_big.m` | Mosek arm, size ladder, and the rungs Mosek cannot solve |
| `bl_fix.m`, `bl_bisdump.m` | pin γ to pose the feasibility question; write the bisection programs |
| `dump2cuadmm.m`, `cuimport.m`, `bl_verify.m` | SDP → cuADMM input; cuADMM certificate → scored against the ORIGINAL data |
| `bl_expect.tsv` | banked expectations — a package input, kept in the repo |
| `gpu/` | the cuADMM arm: `bl_gpu.sh`, `bl_bisgpu.sh`, `bl_gpu_post.sh` |
| `private/` | plant builders (`bl_b_*`) and helpers (`sdpshape`, `stab_mirror`, `opnorm_pi`, …) |
| `experiments/` | the one-off studies behind specific claims in the report; see its README |
| `shadows/` | instrumented copies of two PIETOOLS functions, stored as `.txt`; see below |
| `results/2026-09-24/` | the measured tables and raw logs behind `BASELINE_REPORT.md` |
| `cuadmm_outdir.m`, `cuadmm_private.m`, `cuadmm_shadow.m` | output root; private-helper hook; shadow control |

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

## Shadows: read this before running an experiment that uses one

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
the low-rank ladder, and the psatz claim that has not been isolated for a searched Lyapunov
operator. It also records an **inference from the code, not a measurement**: in 2-D both
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
