# cuADMM testing regime, 2026-09-26

Run on the workstation from **04:56 to 07:44** (2 h 48 min of a window ending 16:15), on fresh
SDP dumps of the current tree. All numbers below were measured in this run unless marked
*inferred*. Driver: `gpu/regime.sh`; blocks: `bl_regime.m`; the bisection: `bl_bisect.m`.
Code: HEAD `5069896c` at start; `2478d265` from B2 on (no lean from the start-up transient);
`b9dac69c` for the second pass (objective-form γ from the dump; B1e). `env.txt` has the dirty
files, which were all outside the executive path (sopvar, lowrank settings). Those three
hashes (and `env.txt`'s) were later combined into the single commit "cuadmm: certified
bisection (bl_bisect), regime harness, 2026-09-26 regime results". The code differences
between the stages are the two changes just named; the per-stage commits remain only in the
local reflog (old tip `752cb47f`) until it expires.

## Verdict

- **cuADMM certified bisection is sound everywhere it was tested.** Every cuADMM certified
  bound is at or above Mosek's certified-infeasible end. The 4 held-out controls at Mosek's
  certified-infeasible γ all read not-F. The 1.01λ\* stability sentinel is correctly not certified.
- **It certifies where cuADMM reaches pinf ≈ 1e-7, and nowhere else.**
  - That holds for 1-D primal H∞ (bounds 1.0002× Mosek's) and for 8 of 9 1-D feasibility
    programs.
  - It fails where cuADMM plateaus at 1e-6..1e-5:
    - coercive H∞ (`hinfco_rd1`), even at 100k iterations;
    - the nonlinear Fisher case, m = 5,093;
    - 2-D stability with the linear psatz generators, m = 4,492.
  - These failures are **accuracy limits, not time limits.**
- **At the sizes tested, Mosek is faster by 2–3 orders of magnitude.** A 1-D certified
  bisection takes Mosek ~1 s and cuADMM 5–15 min.
- **The objective-class ladder puts the time crossover near m ≈ 45k** (*inferred*, from 3 rungs
  of a replicate family). Mosek's own memory wall on 64 GB sits at about the same m.
- **The adoption case is therefore memory and size only.** It rests on reaching the
  accuracy needed to certify at those sizes, which nothing here yet shows.
- **Bisection alone repairs objective-form failures.** Mosek/SeDuMi certify `ctrl_rd1` in
  [128, 1266] (the objective form gave 15433).

## 0. Instrument (B0, B0g, B6b)

- **Dumps.** All 19 fresh dumps reproduce the 09-24 programs exactly (same m; same Mosek rel_b,
  e.g. `hinf_rd1` 5.10481e-07).
- **Binary.** The GPU kernels are in `libcuadmm_lib.so` and `libpsd_lib.so`, both sm_89
  cubins (cuobjdump), built before the 09-24 baseline logs. So the baseline ran on native
  sm_89. The `.sm52.bak` exe links the same new libraries, so it is not an old build.
- **Clean re-times.** Iteration counts are identical to the baseline, and the 09-24 wall
  times were inflated 1.8–3× (overlapping GPU jobs):

  | run | iterations | clean | banked 09-24 |
  |---|---|---|---|
  | `stab_rd1_hv` 1e-6 | 4023 | **5.3 s** | 15.8 s |
  | `stab2_rd` 1e-4 | 1316 | **14.9 s** (11.4 ms/it) | 27.5 s (21 ms/it) |
  | `nl_fisher_nobnd` 1e-4 | 4567 | **72.9 s** (16 ms/it) | 169 s |

- **Objective-form runs.** At 1e-4 with a 120-s cap, none of est/ctrl/hinfco/hinfduco/hinfdu/nl
  converged (the known gap limit). Their clean rates: 2.2 ms/it for est/hinfco/hinfdu,
  1.9 ms/it for ctrl, 0.53 ms/it for hinfduco and 18 ms/it for nl.

## 1. Reference brackets, Mosek and SeDuMi (B1a)

Certified [γ_I, γ_F] from `bl_bisect` (Mosek at rtol 1e-6, SeDuMi at 1e-4). The whole block took 1 min.

| case | Mosek | SeDuMi | objective form (09-24) |
|---|---|---|---|
| `hinf_rd1` | [0.182545, **0.182626**] | [0.182042, 0.183170] | 0.182626 |
| `hinf_rd1_hv` | [0.182401, **0.182482**] | [0.181434, 0.182577] | 0.182482 |
| `hinfdu_rd1` | [0.182470, **0.182556**] | [0.182039, 0.182652] | 0.182551 |
| `hinfco_rd1` | [0.186497, **0.192305**] | [0.189261, 0.191377] | 0.193002 (numerr 2) |
| `hinfduco_rd1` | [76.8, none ≤ 1229] | [38.4, none] | 8251.9 (numerr 2) |
| `ctrl_rd1` | [128, **3917**] | [32, **1266**] | 15433 (numerr 2) |
| `est_rd1` | [—, 4.10e-4] | [—, 4.77e-4] | 2.11e-4 (UNKNOWN) |
| `est_rd1`, pin row rescaled | [—, **3.19e-4**] | [—, 5.08e-4] | |

- **Bisection turns three objective-form failures into certified statements.**
  - `ctrl_rd1`: combining both solvers, γ\* ∈ [128, 1266].
  - `hinfduco_rd1`: certified infeasible up to 76.8 (Mosek). So the 8251.9 may be a genuine,
    very conservative bound rather than nonsense.
  - `est_rd1`: certified feasible at 3.19e-4, but nothing certifies infeasibility anywhere.
- **Near γ\* Mosek often returns UNKNOWN rather than a certificate.** So the certified lower
  end stops short (0.99956 γ\* on `hinf_rd1`) even at rtol 1e-6.
- **SeDuMi's upper ends are looser, but it certifies `ctrl_rd1` 3× tighter than Mosek.**

## 2. cuADMM certified bisection vs Mosek (B1c, B1d, B1e; acceptance run)

Settings: cert_tol 1e-7, psd_tol 1e-7 (held fixed after the hinf_rd1 calibration), psd_abs 1e-6.

| case | cuADMM certified γ_F | / Mosek γ_F | solver time | iterations | held-out at Mosek γ_I |
|---|---|---|---|---|---|
| `hinf_rd1` (acceptance, 09-26) | 0.182663 | **1.00020** | 552 s | 245k | — |
| `hinf_rd1_hv` | 0.182519 | **1.00020** | 660 s | 246k | not F ✓ |
| `hinfdu_rd1`, caps 20k | 0.185686 | 1.0172 | 756 s | 333k | not F ✓ |
| `hinfdu_rd1`, caps 60k (B1e) | 0.182846 | **1.00159** | 911 s | 406k | |
| `hinfco_rd1` | none (5 seeds to 20 γ\*) | — | 319 s | 150k | not F ✓ |
| `hinfco_rd1`, 100k it at 1e-8 (B1e) | none | — | 223 s | 100k | |
| `ctrl_rd1` (B1d, before the fix) | 4700 (seed only) | 1.20 | 1125 s | 661k | not F ✓ |
| `ctrl_rd1`, pin rescaled (B1e) | 57562 | 14.7 | 210 s | 150k | |

- **Soundness.** Every cuADMM F sits at or above Mosek's certified I, and every held-out
  control is not F. cuADMM produced **no exact Farkas certificate** anywhere, so its certified
  lower end is always the analytic or Mosek one.
- **`hinfdu_rd1`.** The dual form converges slowly near γ\*. At 20k caps, feasible probes
  leaned i, and the bound came out 1.7% loose. 60k caps close this to 0.16% for +20% time.
- **`hinfco_rd1`.** The iterate plateaus at pinf ≈ 7.5e-6 even at 100k iterations. The
  repaired λ_min is −1.66e-7 relative and −6.8e-6 absolute, so the point cannot be
  certified.
- **`ctrl_rd1`, two failures of steering, neither of validity:**
  - (B1d) ten probes "converged" at iterations 6–7 inside the start-up transient, with pinf
    2.8e-2 < τ = 3.5e-2. Five of those were below Mosek's certified-infeasible 128. **Fixed**:
    leans now need ≥ 200 iterations.
  - (B1e) with `pin_scale 'auto'`, the pin row's coefficient is ‖b‖/γ_hi ≈ 3e-4. cuADMM
    does not scale small rows up, so the pin barely constrains the solver. Seeds from 4700 to
    37603 never certified in 30k iterations.
  - **Rescaling the pin row helps when γ ≪ ‖b‖ (`est_rd1`) and hurts when γ ≫ ‖b‖.** The fix
    that works both ways is to eliminate the pinned variable instead of appending a row
    (not implemented).
- **Near-exact Farkas vectors (open decision).** On `ctrl_rd1`, cuADMM's y fails the exact
  NSD test only by 5e-15 relative, at b'ŷ ≈ 0.04. That excludes every feasible X with
  trace below ~1e12. Accepting such "bounded-radius" certificates would give cuADMM a
  certified lower end there. It changes the definition of I, so it is left for the
  maintainer.

## 3. Feasibility classes (B2) and stability sentinels (B2S)

One certification run per program (cuADMM at tol 1e-7 with a 50k cap; Mosek), verified by repair.

| case | Mosek | cuADMM | cuADMM time | iterations | cuADMM repaired λ_min (rel / abs) |
|---|---|---|---|---|---|
| `stab_rd1` | F | **F** | 15.6 s | 19,111 | −1.6e-8 / −7.7e-8 |
| `stabdual_rd1` | F | **F** | 15.6 s | 19,111 | same program as `stab_rd1` |
| `stabpde_rd1` | F | **F** | 36.6 s | 50,000 | −4.5e-8 / −1.5e-7 |
| `stabpded_rd1` | F | **F** | 36.6 s | 50,000 | same as `stabpde_rd1` |
| `stab_rd1_hv` | F | **F** | 19.2 s | 14,092 | −8.8e-9 / −3.1e-8 |
| `stab_rd1_tight` (0.95λ\*) | F | **F** | 51.9 s | 37,456 | −1.5e-9 / −2.1e-8 |
| `stab_tr1` | F | **F** | 1.5 s | 1,973 | −5.0e-8 / −4.1e-8 |
| `stab_wave1` | F | **F** | 39.8 s | 26,396 | −6.4e-8 / −5.5e-8 |
| `wellposed_rd1` | F | f | 0.8 s | 928 | −1.18e-7 / −5.0e-8 (misses 1e-7) |

Mosek takes 0.4–1 s per case.

- **Sentinels.**
  - The heavy 1-D reaction–diffusion LPI is Mosek-certified feasible at every fraction up to
    **1.000 λ\*** (λ_LPI ≥ 0.99990 λ\*). The LPI's strictness margin is below resolution, so
    no "just past λ_LPI" sentinel exists.
  - At **1.01 λ\*** Mosek returns a verified Farkas I, and cuADMM is correctly **not** F
    (repaired λ_min −1.5e-5 relative, 150× over the tolerance).
  - No soundness trip.

## 4. Nonlinear, PIESOS Fisher at R = 4, m = 5,093 (B3)

- **Objective form** (solved from the dump, since the builder records no γ): γ_obj = 0.99999.
- **Mosek bisection.** Certified F down to **0.999974**, with **no** certified I (UNKNOWN
  below γ\*). 15 probes, 702 s.
- **cuADMM at 1.2 γ_F cannot be certified:**

  | tol | iterations | time | pinf | repaired λ_min rel / abs |
  |---|---|---|---|---|
  | 1e-4 | 4,490 (converged) | 78 s | 9.4e-5 | −3.0e-6 / −3.5e-6 |
  | 1e-5 | 15,000 (cap) | 261 s | 2.4e-5 | −1.05e-6 / −1.3e-6 |
  | 1e-6 | 25,000 (cap) | 437 s | 1.6e-5 | −7.3e-7 / −9.5e-7 |

  For comparison, Mosek's point repairs to −3.3e-9. B3c (a cuADMM bisection) was skipped by
  its gate.

## 5. 2-D stability with the four linear psatz generators, 0.5 λ\*, m = 4,492 (B5)

- **Mosek.** F, verified (repaired λ_min −8.7e-12, residual 4e-14). The solve took 64 s, and
  the item 120 s including the repair.
- **cuADMM.** **29.7 ms/it** (first measurement on linear4). After 10k iterations (295 s)
  pinf is 8.1e-6 and the repaired λ_min is −7.8e-6: not F.

## 6. Objective-class ladder, `scale_hinf` n = 2, 4, 8 (B4)

| n | m | Mosek bisection (rtol 1e-3) | Mosek bracket | cuADMM ms/it |
|---|---|---|---|---|
| 2 | 1,066 | 2.9 s | [0.182155, 0.182525] | 4.6 |
| 4 | 4,258 | 37.2 s | [0.182155, 0.182525] | 6.1 |
| 8 | 17,026 | 525.8 s | [0.180192, 0.182525] | 9.8 |

- **Scaling.** Mosek's bisection time grows ~m^1.9, cuADMM's per-iteration cost ~m^0.3.
- **Crossover** (*inferred*). A cuADMM bisection at `hinf_rd1_hv`'s ~245k iterations costs
  ~40 min at n = 8, against Mosek's 8.8. They would cross near n ≈ 13 (m ≈ 45k), which is
  about Mosek's memory wall on 64 GB.
- **Caveat.** The family is n uncoupled copies, so iteration counts are flat by
  construction.

## 7. Fixes made during the run (each measured, committed)

1. **The relative PSD test alone was unsound** (shakedown, `acceptance/`).
   - Near γ\*, SeDuMi returned points with λ_max ≈ 3e3 (500× a genuine point). Their
     relative λ_min passed at −4e-8 while the absolute value was −1e-4. Mosek's exact Farkas
     certificates contradict those F verdicts.
   - Fix: F also needs λ_min ≥ −1e-6 absolute and ‖X‖ ≤ 50× the seed's certified ‖X‖. After
     the fix, the SeDuMi and Mosek brackets are consistent.
2. **No lean from runs stopped inside the first 200 iterations** (`ctrl_rd1`, above).
3. **Objective-form γ computed from the dump** when the builder records none (nonlinear).

## 8. Next steps (not done)

- **Eliminate the pinned variable instead of appending a row.** This removes the pin-scaling
  trade-off and the pin-check amplification (§2).
- **Decide on bounded-radius Farkas certificates** for cuADMM (§2). This is the maintainer's
  call.
- **Certification accuracy is the binding limit at m ≳ 5k.** Candidates:
  - facial reduction (the pinned programs have no interior in some Gram blocks);
  - an operator-level margin (eppos) instead of psd_tol;
  - GPU-seed → CPU-refine (memory `gpu-seed-cpu-refine`).
- **Raise kmax near the boundary for slow classes.** `hinfdu_rd1` needed 60k.
- **Correct `BASELINE_REPORT.md`.** Its "cuADMM wins on accuracy on 2-D stability at 0.5 λ\*"
  rests on a residual-grade point (rel_b 2e-4) that this standard would not certify.

## Files

`regime.tsv` (one row per item) · `retime.tsv` · `driver.log` · `env.txt` · `logs/B*.log`
(per block) · `probes/*_probes.tsv` (every probe: verdict, repair diagnostics, iterations,
times) · `acceptance/` (the hinf_rd1 acceptance runs and the SeDuMi guard shakedown). The
per-run X and y (541 MB) stay in `C:\Users\mpeet\cuadmm_out\regime_0926\regime\bisect`.
