# cuADMM harness: Sol-readiness regime, 2026-09-27

Workstation (RTX 4090, native sm_89; binary and library sha256 in `env_desktop.txt`). Output
folder `C:\Users\mpeet\cuadmm_out\regime_0927`. Driver `gpu/regime.sh`, blocks in `bl_regime.m`.

- **Night 1** ran 04:32–06:47 (blocks B0g B1a B2S B2 N2b S1 N4 B1c B1d).
- **S2** (the coupled ladder) ran 06:47–11:13.
- **Shakedown**, before either: B0 and N1 at 03:39.

Everything ran with the default acceptance rule `psd_clip`: η of the repaired point, clipped to
PSD, must be ≤ 1e-7, plus the size guards. Every number below is measured unless it says
*inferred*.

## 1. Soundness: nothing infeasible was accepted

| test | what must happen | result |
|---|---|---|
| sentinels `sent101` (1.01 λ\*), `sentLPI` (1.000 λ\*) | not F | not F (50k it) |
| held-out controls at Mosek's certified γ_I (hinf_rd1_hv, hinfdu, hinfco, ctrl) | not F | not F, all four |
| **N2b**: bisection seeded at 0.99 and 0.999 × Mosek γ_I, no doubling (the first certification, which has no norm reference) | never F | never F: SeDuMi I or f, cuADMM i, on hinf_rd1 and hinfdu_rd1 |
| 1.01 λ\* decoupled sentinels at n = 16/24/32 (m up to 138,240) | not F | not F (30k it) |
| 1.01 λ\* coupled sentinels at n = 16/24/32 | not F | not F (time limit; see §4) |

No `SOUNDNESS_TRIP` and no `DRIVER_TRIP`. N1 driver checks: 4/4 PASS (`N1_checks.txt`).

## 2. Reproducibility against 09-26

- **Fresh dumps (B0):** all 19 have m, nnz and Mosek rel_b bit-identical to 09-26.
- **Mosek brackets (B1a):** identical to 09-26 on all 8 objective cases.
- **cuADMM bisections (B1c) under the new rule:** the same bounds to every digit.
  - `hinf_rd1_hv` γ_F = 0.1825185127, 1.0002× Mosek.
  - `hinfdu_rd1` γ_F = 0.1856860509, 1.017× at 20k caps.
  - `hinfco_rd1`: still no upper end.
- **Feasibility (B2):** the same 8/9 F and the same iteration counts. `wellposed_rd1` is f again.
- **`ctrl_rd1` (B1d):** now 4,333.7 in 9 s, against Mosek's 3,917 (10.6% loose; 09-26: 4,700 in
  1,146 s). With the unscaled pin row, cuADMM converges in about 30 iterations. The 200-iteration
  lean guard then voids the leans, including true-feasible probes at 3,953 and 4,310. Sound, but
  loose. It is on the problematic list, and the fix is the pinned-variable decision.

## 3. Phase costs for sizing Sol jobs

| case | m | iterations | solve | ms/it | setup | certify | MATLAB memory | verdict |
|---|---|---|---|---|---|---|---|---|
| `scale_stab_n08` | 8,640 | 14,092 | 120 s | 8.5 | 0.1 s | 0.4 s | 2.1 GB | F |
| `scale_stab_n16` | 34,560 | 14,092 | 149 s | 10.5 | 0.2 s | 0.7 s | 2.2 GB | F |
| `scale_stab_n24` | 77,760 | 14,092 | 242 s | 17.2 | 0.2 s | 2.1 s | 2.3 GB | F |
| `scale_stab_n32` | 138,240 | 14,092 | 355 s | 25.2 | 0.4 s | 4.2 s | 2.3 GB | F |
| decoupled sentinel n32 | 138,240 | 30,000 (cap) | 754 s | 25.1 | 0.4 s | 51 s | 2.6 GB | i |
| `Th_n3`, live, tol 1e-4 | 80,550 | 3,063 | 1,058 s | 345 | 238 s | 529 s | 3.4 GB | f (η 5.4e-5) |
| `Th_n3`, kept 1e-6 point, re-certified | 80,550 | (72k, 09-17) | — | — | — | 506 s | — | f (η 2.8e-7) |
| coupled `scale_rot_f0p50_n08` | 8,640 | 21,984 | 226 s | 10.3 | 1.0 s | 0.5 s | — | F (η 2.8e-8) |

- **Decoupled ladder:** 14,092 iterations at every n, the same as n = 1. It tests cost per
  iteration, not convergence: ms/it grows 3× over 16× in m, and the host phases are negligible.
- **Th_n3 is host-bound.** Setup plus certification is about 13 minutes per probe, and on Sol the
  A100 is billed through it. A 1e-4 probe's certification is predictably futile (η 5e-5). The
  kept 1e-6 point from 6.8 h of solving misses the tolerance by 2.8×.
  - *Proposed, not implemented:* skip certification when pinf_min > 100 × psd_eta_tol. Measured
    clipped η has always been ≥ about 0.3 × pinf, so the skip cannot lose a certificate.

## 4. The coupled ladder: cuADMM's cost is the fill of chol(AAᵀ), not m

`scale_rot` (built by `bl_big(ns,frac,true)`) is `x_t = x_ss + A x` with `A = Q diag(μ) Qᵀ` and Q
dense orthogonal. Its answer follows from an orthogonal-invariance argument, and **Mosek confirms
it at n = 2, 4, 8 and 16** (n16: 66 s). cuADMM results:

| n | m | nnz(At) | cuADMM | ms/it |
|---|---|---|---|---|
| 8 | 8,640 | 204k | F, 21,984 it | 10.3 |
| 16 | 34,560 | 1.35M | 12,900 it in the 1,420 s limit, pinf 2.4e-7, not finished | ~110 |
| 24 | 77,760 | 4.24M | 2,700 it in 2,465 s | ~380 |
| 32 | 138,240 | 9.66M | no iteration printed in 3,772 s | — |

The reason is measured symbolically (`aat_fill.tsv`, AMD ordering). The factor of AAᵀ, which
cuADMM builds once and solves against every iteration:

| n | factor nnz, decoupled | factor nnz, coupled | coupled, share of dense |
|---|---|---|---|
| 8 | 0.22M | 3.6M | 9.5% |
| 16 | 0.91M | 54M | 9.1% |
| 24 | 2.1M | 260M | 8.6% |
| 32 | — | 812M (~10 GB) | 8.5% |

Nested dissection is 15–19% *worse* than AMD, so building CHOLMOD with METIS would not help.
Time per iteration tracks the factor size.

**Consequences:**
- The 09-26 "crossover near m ≈ 45k" was inferred from decoupled families, which keep AAᵀ
  block-diagonal. **It does not apply to coupled problems.**
- On coupled problems cuADMM still needs far less memory than an IPM (about 10 GB at n = 32
  against Mosek's 142 GB dense Schur complement). But its per-iteration cost grows like the
  factor, about m², and the one-time factorization and triangular-solve setup at n = 32 did not
  finish within an hour here.
- *Inferred Sol costs, A100 at 1.45–2× the 4090:*
  - coupled n = 16 to a certificate: 20–25k iterations, 25–35 min, **15–20 CHE**;
  - coupled n = 24: about 2 h, **~70 CHE**;
  - coupled n = 32: 5+ h once past setup, **> 180 CHE**.

  Each needs your approval.
- *Open question:* can cuADMM solve the AAᵀ system iteratively (preconditioned CG) instead of by
  Cholesky? That is what would change the coupled-problem scaling.

## 5. Sol-0

Prepared (`sol/sol0.sh`), not run: the workstation is not on the ASU VPN. The desktop reference
for the same blocks and dumps (`x0_desktop_*`):
- N1: 4/4 PASS;
- `stab_tr1`: 1,973 it, 1.57 s, F;
- `stab_rd1_hv`: 14,092 it, 19.0 s, F;
- `scale_rot_f0p50_n04`: 10,894 it, 52.2 s, F.

## Files

- `regime.tsv`, `driver.log`: manifest and driver log.
- `probes/`: every run's probe rows, including t_init, t_cert, mem_mb and F_etapsd.
- `bl_mosek.tsv`: fresh Mosek references with the η columns.
- `bl_big.tsv`: builds without solving.
- `calpsd.tsv`: the 120-iterate calibration of psd_eta_tol.
- `aat_fill.tsv`: the factor-fill measurements.
- `x0_desktop_*`: the Sol-0 desktop reference.
