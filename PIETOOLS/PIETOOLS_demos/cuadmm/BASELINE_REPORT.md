# PIETOOLS baseline: Mosek vs cuADMM across all application classes

Measured 2026-09-24 on one workstation. 37 cases spanning every executive class PIETOOLS
ships, plus a size ladder that runs past what an interior-point method can do here. The
purpose is a reference point per **application class** — not per spatial dimension — so a GPU
first-order solver can be judged on a tradeoff curve rather than on anecdotes.

## Method

Each case builds a self-contained plant inline from the `pde_var` API (never the
Examples_Library, whose files run `evalin('base',...)` side effects and at least one of which
prompts with `input()`), calls the stock executive, and records the SDP it actually solved.

- **Capture is post-solve.** Every executive calls `lpisolve` internally, and for programs
  carrying an `ineq` expression `sossolve` inserts blocks via `addextrasosvar` *during* the
  solve (`sossolve.m:177-183`). A pre-solve capture of an H-infinity program misses the
  `gamma>=0` block.
- **`sol_opts.simplify` is forced false.** `sossolve.m:256-262` overwrites `At/b/K/c` with the
  `sospsimplify`-reduced system, so with it on Mosek would not be solving the SDP that gets
  dumped. Verified a no-op on the PIESOS cases; not assumed elsewhere.
- **The dump is validated.** `rel_b` (Mosek inside `sossolve`) is compared with `rel_dump`
  (Mosek re-solving the dumped bytes through `Sedumi2Mosek`). Of the 35 successful rows, 31
  agree to every printed digit, 3 agree to better than 1e-4 relative, and one does not
  reproduce: `hinf2du_rd`, 4.84e-03 against 5.17e-03 (6.5%). That row has numerr = 2, and
  Mosek's answer there is not stable from one solve to the next. Mosek's own optimizer clock is
  `solinfo.info.cpusec`, from `MSK_DINF_OPTIMIZER_TIME`.
- **cuADMM results are verified, not trusted.** `cuimport.m` rebuilds each certificate in
  SeDuMi cone coordinates and re-scores it against the pre-dump `At/b/c`. Across 34
  certificates the rebuilt Gram blocks are exactly symmetric (`asym = 0`). The verified row
  residual is ~2x **cuADMM's own reported primal infeasibility**, which is precisely the
  documented `(1+||b||)` vs `||b||` denominator difference, and so at most ~2x the requested
  tolerance. Some rows sit well below the tolerance: `hinf_rd1_hv` reaches 3.65e-05 at 1e-4,
  and `stab_wave1` reaches 3.86e-07 at 1e-6. On these problems the solver's self-report is
  honest.
- **cuADMM convergence is judged by the presence of a `CUADMM_TIMING` line**, never by exit
  code: `cuadmm_exe` exits 0 on total CUDA failure, and a signal-killed run exits 143.

## Hardware, and three caveats that bound every number

Intel i9-14900KF (24C/32T), 64 GB RAM, RTX 4090 24 GB. Mosek gets 24 cores and 64 GB; cuADMM
gets one GPU. That asymmetry is the honest comparison — each solver on the hardware it targets
— but no ratio should be quoted without it.

1. **The workstation was shared.** Another session's MATLAB was resident for part of the run.
   `nl_fisher_opt` measured 48.3 s and 33.5 s on two passes of the *identical* program, so
   treat individual times as +/-30%.
2. **Part of the GPU arm ran twice concurrently.** Killing the first orchestrator killed its
   bash wrapper but not its `wsl.exe` child, so two `cuadmm_exe` instances shared the GPU for
   the alphabetical range `est_rd1`..`nl_fisher_*`. **Objective-case GPU times in that range
   are inflated by an unknown factor and should be read as upper bounds.** The scaling ladder,
   the `stab_*` rows and the `wellposed` row are single-entry and clean. The table below was
   rebuilt from the per-case logs rather than from the appended results file, so each
   (case, tolerance) appears once.
3. **cuADMM runs were wall-capped** (120 s at 1e-4, 300 s at 1e-6 in the main sweep; larger
   caps on the re-runs). `no` below means "did not converge within the cap given", which is
   not the same as "cannot converge" — `hinf_rd1` shows `no` at 1e-4 under a 120 s cap yet
   converged at 1e-6 given 220 s.

## Master table

`m` constraints, `K.f` free variables, `norm_b` is `||b||` after assembly, Mosek time is its
own optimizer clock, `rel_b = ||At'x-b||/||b||`. cuADMM columns are solve seconds, `no` (did
not converge in the cap) or `—` (not attempted). `cu_gap` is the relative duality gap reached
at 1e-4 — the quantity that decides the whole comparison. `verif` is the independently
recomputed row residual of cuADMM's own 1e-4 certificate.

| case | class | kind | m | K.f | blocks | norm_b | mosek s | rel_b | gam | numerr | cu 1e-4 | cu 1e-6 | cu_gap | verif |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| stab_rd1 | stability | feas | 130 | 43 | [10 15 6] | 9.43e-07 | 0.0104 | 1.32e-08 | — | 0 | 1.2 | 16.0 | 6.74e-07 | 0.000194 |
| stabdual_rd1 | stability | feas | 130 | 43 | [10 15 6] | 9.43e-07 | 0.00819 | 1.33e-08 | — | 0 | 1.2 | 16.0 | 6.74e-07 | 0.000194 |
| stabpde_rd1 | stability | feas | 54 | 0 | [10 15 6] | 8.03e-06 | 0.00899 | 5.09e-08 | — | 0 | 6.5 | 111.8 | 0.0001 | 6.5e-05 |
| stabpded_rd1 | stability | feas | 54 | 0 | [10 15 6] | 8.03e-06 | 0.00941 | 5.1e-08 | — | 0 | 9.5 | 132.9 | 0.0001 | 6.5e-05 |
| stab_rd1_hv | stability | feas | 135 | 45 | [11 20 11] | 9.43e-07 | 0.00781 | 2.79e-08 | — | 0 | 2.0 | 15.8 | 4.95e-06 | 8.96e-05 |
| stab_rd1_tight | stability | feas | 135 | 45 | [11 20 11] | 9.43e-07 | 0.00928 | 1.15e-08 | — | 0 | 4.6 | 31.3 | 9.82e-06 | 0.000198 |
| stab_tr1 | stability | feas | 96 | 43 | [10 15 6] | 2e-06 | 0.00435 | 6.42e-08 | — | 0 | 0.5 | 2.0 | 3.63e-05 | 0.000155 |
| stab_wave1 | stability | feas | 520 | 172 | [20 30 12] | 1.33e-06 | 0.021 | 5.9e-08 | — | 0 | 1.5 | 37.2 | 7.18e-09 | 0.000149 |
| hinf_rd1 | hinf | obj | 242 | 44 | [10 31 13 1] | 0.707 | 0.0222 | 5.1e-07 | 0.183 | 0 | no | 219.7 | 0.000259 | — |
| hinfco_rd1 | hinf | obj | 176 | 1 | [10 31 13 1] | 0.707 | 0.116 | 9.39e-07 | 0.193 | 2 | no | no | 0.00237 | — |
| hinfdu_rd1 | hinf | obj | 244 | 44 | [10 31 13 1] | 1.41 | 0.0344 | 1.52e-08 | 0.183 | 0 | no | no | 0.00343 | — |
| hinfduco_rd1 | hinf | obj | 168 | 1 | [10 31 13 1] | 1.41 | 0.0461 | 0.00186 | 8.25e+03 | 2 | no | no | 0.376 | — |
| hinf_rd1_hv | hinf | obj | 267 | 46 | [11 41 23 1] | 0.707 | 0.0257 | 4.83e-08 | 0.182 | 0 | 72.4 | no | 9.98e-05 | 3.65e-05 |
| h2c_rd1 | h2 | obj | 272 | 44 | [10 1 30 16 12 7 1 1] | 1.41 | 0.0327 | 2.78e-07 | 0.289 | 0 | no | no | 0.00173 | — |
| h2cco_rd1 | h2 | obj | 56 | 1 | [10 15 6 1 1] | 1 | 0.0262 | 1 | 1.01 | 0 | no | no | 1 | — |
| h2o_rd1 | h2 | obj | 269 | 44 | [10 8 16 30 7 12 1] | 0.707 | 0.0626 | 3.97e-06 | 0.00263 | 2 | no | no | 0.0171 | — |
| h2oco_rd1 | h2 | obj | 56 | 1 | [10 15 6 1 1] | 0.5 | 0.0098 | 5.68e-08 | 9.03e-05 | 0 | no | no | 0.000456 | — |
| est_rd1 | estimator | obj | 189 | 4 | [10 32 14 1] | 0.707 | 0.0648 | 1.07e-06 | 0.000211 | 0 | no | no | 0.00528 | — |
| h2est_rd1 | estimator | obj | 198 | 25 | [10 17 30 8 12 1 1] | 0.707 | 0.0696 | 6.45e-07 | 0.000171 | 2 | no | no | 0.0152 | — |
| ctrl_rd1 | controller | obj | 181 | 4 | [10 32 14 1] | 1.41 | 0.0494 | 0.00108 | 1.54e+04 | 2 | no | no | 0.432 | — |
| h2ctrl_rd1 | controller | obj | 202 | 8 | [10 30 17 12 8 1 1] | 1.41 | 0.0476 | 0.000781 | 1.84e+04 | 2 | no | no | 0.406 | — |
| wellposed_rd1 | wellposed | feas | 108 | 0 | [10 10 15 15 6 6] | 17.7 | 0.00983 | 8.74e-08 | — | 0 | 0.2 | 1.2 | 9.14e-08 | 0.000174 |
| nl_fisher_opt | nonlinear | obj | 5093 | 1 | [5 50 50 40 456 380] | 18.5 | 33.5 | 5.11e-06 | — | 0 | no | no | 0.239 | — |
| nl_fisher_gamfix | nonlinear | feas | 5093 | 0 | [5 50 50 40 456 380] | 18.5 | 36.4 | 5.15e-06 | — | 0 | 186.1 | no | 6.64e-05 | 0.0002 |
| nl_fisher_nobnd | nonlinear | feas | 4978 | 0 | [5 50 456 380] | 18.4 | 39.1 | 3.44e-05 | — | 0 | 168.9 | no | 6.88e-05 | 0.000192 |
| stab2_rd | stability | feas | 3456 | 0 | [8 424] | 1.47e-05 | 13.6 | 2.08 | — | 2 | 27.5 | no | 5.34e-05 | 0.0002 |
| stab2dual_rd | stability | feas | 3456 | 0 | [8 424] | 1.47e-05 | 16.7 | 0.00843 | — | 2 | 27.6 | no | 5.34e-05 | 0.0002 |
| hinf2_rd | hinf | obj | 14005 | 1 | [8 869] | 0.5 | 283 | 0.000136 | 1.27 | 2 | no | no | 0.0502 | — |
| hinf2nc_rd | hinf | obj | 16949 | 42 | [8 1229] | 0.5 | 574 | 0.000154 | 3.93 | 2 | no | no | 0.0313 | — |
| hinf2du_rd | hinf | obj | 14006 | 1 | [8 869] | 1.41 | 307 | 0.00484 | Inf | 2 | no | no | 0.355 | — |
| h2_2dc_rd | h2 | obj | 4127 | 1 | [8 614 1] | 1 | 46.3 | 1 | 1.1 | 0 | no | no | 1 | — |
| scale_stab_n02 | scaling | feas | 540 | 180 | [22 40 22] | 1.33e-06 | 0.0373 | 2.33e-08 | — | 0 | 5.9 | 45.4 | 4.95e-06 | 8.96e-05 |
| scale_stab_n04 | scaling | feas | 2160 | 720 | [44 80 44] | 1.89e-06 | 0.2 | 3.87e-08 | — | 0 | 9.1 | 70.0 | 4.95e-06 | 8.96e-05 |
| scale_stab_n08 | scaling | feas | 8640 | 2880 | [88 160 88] | 2.67e-06 | 3.82 | 6.29e-08 | — | 0 | 7.9 | 59.7 | 4.95e-06 | 8.96e-05 |
| scale_stab_n16 | scaling | feas | 34560 | 11520 | [176 320 176] | 3.77e-06 | 118 | 5.84e-08 | — | 0 | 12.6 | 82.5 | 4.95e-06 | 8.96e-05 |
| scale_stab_n24 | scaling | feas | 77760 | 25920 | [264 480 264] | — | **n/a** | — | — | — | 19.9 | — | 4.95e-06 | 8.96e-05 |
| scale_stab_n32 | scaling | feas | 138240 | 46080 | [352 640 352] | — | **n/a** | — | — | — | 28.4 | — | 4.95e-06 | 8.96e-05 |

The last two rows are the ladder rungs Mosek cannot solve on this machine: their SDPs
were built and dumped, but only cuADMM produced a certificate. `scale_*` uses heavy
settings, everything else light unless the id says `_hv`; `gam` is blank for
feasibility cases by construction.

## 1. The finding that decides cuADMM's usefulness here

**Viability splits on objective-form vs feasibility-form, and the split is the duality gap.**

Every one of the 18 objective-form cases is gap-limited, and most reach feasibility long
before they run out of time:

| case | m | final pinf | final dinf | final relgap |
|---|---|---|---|---|
| `hinf_rd1` | 242 | 1.95e-05 | 1.52e-05 | **2.59e-04** |
| `h2oco_rd1` | 56 | 1.36e-05 | 1.42e-05 | **4.56e-04** |
| `hinf_rd1_hv` | 267 | 1.83e-05 | 1.69e-05 | 9.98e-05 (just converged) |

pinf and dinf are already inside 1e-4; the gap is not, and on `h2oco_rd1` it falls by about
1e-8 per 100 iterations. Meanwhile **all 19 feasibility-form cases converged**, limited by
pinf or dinf with the gap several orders below tolerance.

Only 2 of 18 objective cases converged (1 of 18 at 1e-4) at any tolerance (`hinf_rd1_hv` at 1e-4 in 72 s;
`hinf_rd1` at 1e-6 in 220 s and 39,122 iterations). Objective cases need roughly 10-40x more
iterations than feasibility cases of comparable size.

**Practical consequence.** Bisection on gamma turns an objective LPI into a sequence of
feasibility LPIs — exactly the form cuADMM finishes. It is already the preferred method for
accuracy reasons; it is also the only form in which this solver is usable.

A second, smaller result points the same way: the *heavy*-settings `hinf_rd1_hv` converged
where its light twin `hinf_rd1` did not, at the same tolerance and a larger SDP. Heavier
settings converging faster is now measured twice.

## 2. The scaling ladder, and where the crossover is

1-D stability, n decoupled states, heavy settings. `m = 135 n^2` exactly. These rows are
clean (no concurrency contamination).

| n | m | nvar | Mosek | cuADMM 1e-4 | GPU mem |
|---|---|---|---|---|---|
| 2 | 540 | 2,748 | 0.037 s | 5.90 s | 806 MiB |
| 4 | 2,160 | 10,992 | 0.200 s | 9.13 s | 808 MiB |
| 8 | 8,640 | 43,968 | 3.82 s | 7.93 s | 814 MiB |
| 16 | 34,560 | 175,872 | **117.8 s** | **12.62 s** | 850 MiB |
| 24 | 77,760 | 395,712 | *45 GB dense Schur — not solvable here* | **19.92 s** | 906 MiB |
| 32 | 138,240 | 703,488 | *142 GB dense Schur — not solvable here* | **28.43 s** | 976 MiB |

- **Crossover between m = 8,640 and m = 34,560.**
- Fitted exponents: Mosek ~1.94 and steepening; cuADMM **0.28** over a 256x range in m.
- **The memory ratio is the real argument, not speed**: under 1 GB of GPU memory at
  m = 138,240 where Mosek's dense Schur complement alone needs 142 GB. About 150x.
- Mosek's ladder ends by *memory*, not patience: an n=24 run reached 22.6 GB resident with
  18.6 GB free and was stopped before it could thrash.
- Iteration counts are identical across all rungs (519 at 1e-4, 4023 at 1e-6). That is the
  documented separability artifact of n decoupled copies, and it doubles as an instrument
  check — as does the verified `rel_b` being identical (8.96e-05) on every rung.

The m = 138,240 certificate, which Mosek cannot produce on this machine, verifies
independently at `rel_b = 9.0e-05` with Gram blocks PSD to 5e-11.

## 3. Defects the baseline exposed (none of them about cuADMM)

**A build wall inside PIETOOLS, now fixed and committed (`0c2bc4b4`).** n=32 could not be
*built*: `lpi_eq` requested a 9,613,344 x 1536 dense array (123.8 GB) at
`if ~all(all(C.C==0))`, because `C.C` is sparse and `sparse == 0` returns a *full* logical of
the same size. Replaced with `nnz(C.C)~=0` — exactly equivalent, O(1). n=32 now builds in
93 s. **Regression:** `bl_check('structure')` passes 22/22 against the pre-fix baseline, with
m, K.f, the block list and nnz(At) all identical (`results/2026-09-24/check_structure.tsv`).
An earlier five-case check (`logs/regress_lpieq.log`) also found identical structure, but its
own verdict reads `false`: `stab_rd1`'s rel_b moved by 2.7e-06 relative, just above that
script's 1e-6 gate. That's Mosek's run-to-run variation, not a change in the program, which is
why `bl_check` uses a 10x band on rel_b instead.

*Not patched, because not measured failing*: the same pattern at
`opvar/2D/@dopvar2d/set.m:472,506` and `dopvar2d.m:588,601`.

**Two programs both solvers reject.** `h2cco_rd1` (H2_norm_c_coercive) and `h2_2dc_rd`
(H2_norm_2D_c) give Mosek `rel_b = 1.000` — and `h2cco_rd1` does it with **numerr = 0**, so
the status word does not catch it — while cuADMM leaves pinf frozen at its *initial* 5.00e-01
with the dual objective running to 2e+10. Two independent solvers failing identically points
at the programs, not the solvers.

**Objective values that are not usable.** `hinfduco_rd1` returns gam = 8251.9;
`ctrl_rd1` 15433; `h2ctrl_rd1` 18396; `hinf2du_rd` returns `Inf`. All carry numerr = 2. The
controller cases also show `normx ~ 7e7`, which is where the post-solve operator inversion in
those executives would bite.

**The one case with a ground truth is 7.4x conservative.** `hinf2_rd` copies the shipped
`Ex_2D_ReactionDiffusion_DDDD` plant, whose example file carries a closed-form L2 gain of
0.1711. The executive returns 1.272 — a valid upper bound, but far from tight, and Mosek
reports UNKNOWN. The other three 2-D objective cases on the same plant return 3.93 and `Inf`.

**On 2-D stability it is MOSEK that fails, not the LPI.** `stab2_rd` (m=3456): Mosek returns
rel_b = 2.08 with numerr = 2 and `||x|| ~ 0.028` — it never satisfied the equalities at all.
cuADMM converges on the same dumped bytes, and its certificate VERIFIES independently at
rel_b = 2.0e-04, with Gram blocks PSD to 3e-10 and exactly symmetric. The dual form tells the
same story (Mosek 8.4e-03, cuADMM 2.0e-04). So the 2-D stability LPI is feasible at
0.5 lam\* and the interior-point solver is what falls over — which corrects an earlier reading
of this row as an LPI failure.

This is the one place in the suite where the GPU solver wins on **accuracy** rather than size,
and it happens at m = 3456, well *below* the crossover. Worth knowing before dismissing a
first-order solver at small m.

## 4. What ||b|| says about the classes

Feasibility cases sit at `||b|| ~ 1e-6..8e-6` — that *is* `eppos2` — and objective cases at
0.5..1.41. Two exceptions matter: `wellposed_rd1` at **17.68** and the nonlinear PIESOS cases
at **18.4**. Well-posedness is the only 1-D feasibility class with an O(1) right-hand side,
which is what separates "the solver handles a tiny b" from "the solver handles this cone". It
is also the fastest case on the GPU by a wide margin: 0.24 s and 44 iterations.

## 5. Nonlinear (PIESOS)

Both feasibility forms converge on the GPU with a realistic cap and PSD violation exactly
zero: `nl_fisher_nobnd` 168.9 s and `nl_fisher_gamfix` 186.1 s, against Mosek's 39.1 s and
36.4 s. Mosek is still 4-5x faster at m ~ 5000, consistent with the crossover sitting above
it. The objective form `nl_fisher_opt` capped, gap-limited, as the split in section 1
predicts.

## 6. Reproducing this

The harness is in this folder (`PIETOOLS_demos/cuadmm`; see `README.md` for the layout):
`bl_cases.m` (registry), `private/bl_b_*.m` (builders), `bl_run.m` (Mosek arm), `bl_scale.m`
and `bl_big.m` (ladder), `gpu/bl_gpu.sh` (GPU arm), `cuimport.m` and `bl_verify.m`
(independent verification), `gpu/bl_gpu_post.sh` (limiter column). The tables behind every
number in this report are in `results/2026-09-24/`, with the raw logs in
`results/2026-09-24/logs/`. The SDP dumps (~1.9 GB) were not kept; re-running the harness
regenerates them under `cuadmm_outdir()`.

This report was measured with the harness in a flat scratchpad layout, before it moved here.
The move rewrote paths only (see `README.md`); it did not change what any case builds or
solves, and the harness has not been re-run in its new location yet.

`bl_run`, `bl_scale`, `bl_verify` and `gpu/bl_gpu.sh` append one row per case and skip ids
already recorded, so they resume where they left off. `bl_run` warns when it does this, because
the resume state persists under `cuadmm_outdir()` and would otherwise mix rows measured on old
code with new ones. `bl_big`, `bl_bisdump`, `gpu/bl_bisgpu.sh`, `gpu/bl_gpu_post.sh` and
`bl_check` start from scratch every time. Two environment notes that cost real time: cuADMM needs
`LD_LIBRARY_PATH=/usr/lib/wsl/lib` (a distro `libcuda.so.595.91.07` shadows the WSL driver
and makes every CUDA call report "no CUDA-capable device"), and `cuadmm_exe` **exits 0 on
total CUDA failure**, so never gate on its exit code.

---

# 7. Bisection as the default for objective classes

Section 1 showed the objective form is gap-limited and that bisection should be the default,
since it turns an objective LPI into a sequence of feasibility LPIs. This section measures
that, on the feasibility question itself.

## The transformation, and why it needs no new code per executive

The stock 1-D executives hard-declare the objective — `dpvar gam; lpidecvar; lpi_ineq;
lpisetobj` (`PIETOOLS_Hinf_gain.m:128-131`) — and only the 2-D ones accept a `gain` argument
(`PIETOOLS_Hinf_gain_2D` even ships `settings.bisect_opts` with min/max/start). The 1-D
estimator's `gain` argument only dispatches to the 2-D version.

Rather than mirror eight executives, `bl_fix.m` pins gamma on the assembled SDP: the objective
vector is a single unit entry, so appending one row `e_j' x = gam` and setting `c := 0` poses
exactly the feasibility question, with no re-derivation to get wrong. b is rescaled *after*
the row is appended, since `sossolve` normalises b.

**Validated against a known answer.** `hinf_rd1` (objective form: gam\* = 0.182626):

| gam/gam\* | Mosek prosta | rel_b | \|\|x\|\| |
|---|---|---|---|
| 0.999 | `PRIMAL_INFEASIBLE` (certificate) | 1.0000 | 3e-18 |
| 1.000 | `PRIMAL_AND_DUAL_FEASIBLE` | 4.6e-07 | 1.4e+02 |

The flip lands on gam\* to six digits. All three 1-D Hinf variants behave identically.

## Per-step cost on the feasibility question

| case | m | Mosek/step | Mosek verdict | cuADMM/step (5000 it) | s/iter |
|---|---|---|---|---|---|
| `h2oco_rd1` | 57 | 0.007-0.010 s | **invalid: pinned the wrong quantity** (see below) | 4.4 s | 0.87 ms |
| `est_rd1` | 190 | 0.038-0.073 s | **UNKNOWN at every gamma** | 10.6 s | 2.13 ms |
| `h2c_rd1` | 273 | 0.018-0.039 s | decides, with certificate | 12.3 s | 2.46 ms |
| `h2est_rd1` | 199 | 0.045-0.091 s | **UNKNOWN at every gamma** | 13.3 s | 2.66 ms |
| `hinf_rd1` | 243 | 0.017-0.029 s | decides, with certificate | — | — |
| `hinf2_rd` | 14006 | **341-480 s** | **UNKNOWN at every gamma** | **172 s** | 34.5 ms |

**The crossover appears on an objective class too.** At m = 14,006 a cuADMM bisection step
(172 s) is 2-2.8x cheaper than a Mosek one (341-480 s) — and Mosek returns UNKNOWN, so the
cost buys nothing there. Consistent with the m ~ 10-30k crossover from the stability ladder.

At small m, Mosek wins by 100-1000x per step, so bisection does NOT rescue cuADMM below the
crossover; it only makes the question well-posed for it.

**Cost of the whole bisection.** For Mosek on the 1-D Hinf class, ~21 steps to 1e-6 relative at
~0.03 s each is about 0.6 s, against 0.022 s for one objective solve — roughly 27x more
expensive. That is the price of a verdict that comes with an infeasibility certificate instead
of a gamma reported alongside numerr = 2.

## The decision rule differs by solver, and cuADMM's is budget-dependent

Mosek returns `PRIMAL_INFEASIBLE_CER` — an actual certificate, no threshold needed.

cuADMM has no infeasibility certificate, and on these problems **neither side converges**,
because the stopping rule needs the duality gap. What separates them is the primal residual
PLATEAU. On `h2c_rd1` the infeasible program sticks at pinf ~ 1.5e-02 from iteration ~500 and
never improves, while the feasible one keeps descending:

| iteration | infeasible pinf | feasible pinf | ratio |
|---|---|---|---|
| 500 | 1.75e-02 | 9.65e-03 | 2x |
| 1000 | 1.62e-02 | 4.55e-03 | 4x |
| 5000 | 1.50e-02 | 5.08e-04 | **30x** |

At a 5000-iteration budget a threshold of 1e-3 separates them cleanly, and **cuADMM's verdicts
then agree with Mosek's certificates on every `h2c_rd1` rung.**

**RETRACTED (2026-09-25): the `h2oco_rd1` rows.** An earlier version of this section said
`h2oco_rd1` was feasible even at 0.80γ\*, with cuADMM's pinf identical to three digits at all
four γ, and concluded that the objective form's γ\* was not the relaxation's true optimum.
That conclusion is wrong, and the cause is mine. `PIETOOLS_H2_norm_o_coercive.m:172` returns
`gam = sqrt(objective)`, so the SDP variable that `bl_fix` pins is γ² ≈ 8.2e-09, not
γ = 9.03e-05. The four rungs pinned it at 0.8–1.2 × 9.03e-05, about 9,000x above its optimum,
so every rung sat deep in the feasible region. That accounts both for "feasible everywhere"
and for the pinf values not depending on γ. The two solvers' agreement there is therefore not
evidence of anything. `bl_bisdump.m` now pins γ² for this case. It has not been re-run.

Only the two H2 *coercive* executives take the square root: the estimators, the controllers and
every H∞ executive return the objective directly (checked). So the `h2c_rd1`, `est_rd1`,
`h2est_rd1`, `hinf_*` and `hinf2_rd` rows are unaffected.

## BUT the budget rule does not transfer to the 2-D case — measured, not assumed

The same rule applied to `hinf2_rd` said every gamma down to 0.18 was feasible, which would
mean bisection recovers a bound 7x tighter than the objective form's 1.272 (the shipped
example's closed-form gain is 0.1711). A control settles it: the Hinf relaxation **cannot** be
feasible below the true gain, so gammas underneath 0.1711 must be rejected.

| gam | gam/exact | pinf @5000 | rule says |
|---|---|---|---|
| 0.05 | 0.29 | 1.18e-03 | infeasible |
| **0.12** | **0.70** | **4.64e-04** | **FEASIBLE — wrong** |
| 0.18 | 1.05 | 2.89e-04 | feasible |
| 0.60 | 3.51 | 5.88e-05 | feasible |
| 1.018 | 5.95 | 2.45e-05 | feasible |

pinf decays smoothly and monotonically straight through the true boundary, with none of the
plateau that made the 1-D case decidable. So at m = 14,006 and 5000 iterations cuADMM is still
in its transient: pinf measures how far it has got, not whether a solution exists. **The
apparent tighter bound is an artefact of the budget and is withdrawn.** Deciding this case
would need a far larger budget — at 34.5 ms/iter, 50,000 iterations is 29 min per step and
~10 h for one bisection.

## What this means for the suite

1. **Use bisection by default on the objective classes.** The objective form is gap-limited
   for cuADMM and returns numerr = 2 (and in four cases a useless gamma: 8251.9, 15433, 18396,
   `Inf`) under Mosek. The feasibility form is well-posed and, where it decides, both solvers
   agree.
2. **Report the per-step cost, and the step count separately.** They scale differently: the
   step count is set by the bisection tolerance, the per-step cost by m.
3. **Two classes are not decidable by either solver as posed** — the Hinf and H2 estimators,
   where Mosek returns UNKNOWN at every gamma. That is a property of those programs, not of
   the method, and it needs diagnosis before either solver is judged on them.
4. **A fixed-iteration threshold is not a portable decision rule.** It was validated at
   m ~ 270 and refuted at m = 14,006 by a control that used a known analytic answer. Any
   cuADMM-driven bisection needs a stagnation test, not a budget — and needs a control of this
   kind on each new problem size before its verdicts are trusted.
