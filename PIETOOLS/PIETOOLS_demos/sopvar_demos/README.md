# heatNd: an N-D heat-equation stability benchmark (N = 1, 2, 3)

Initial coding MMP, 09/27/2026. Revised the same day after two reviews (kappa
parameterization, certified-bracket bisection, reproducible dumps; see the file
headers). §11 (added the same day) is the PIE-form Poincaré diagnostic of the heat LPI,
which runs separately from it. Revised again after the final reviews of both (MMP,
09/27/2026): the Farkas gate is on the absolute violation, a certified +1 above λ1 stops
the bisection, test (9) of the diagnostic is a Rayleigh–Ritz check, and the 3-D
diagnosis states its resolution. No reported bracket changed.

Library face codes (MMP, 09/27/2026). `copquadvar` now takes `options.psatz` = 2d+1 / 2d+2,
the linear weight of one face of the box. The Psatz face terms of `heatNd_lpi` and
`heatNd_poincare` now come from `poscopvar` with these codes. Before, they came from
`heatNd_posw`, a copy of `copquadvar` with a weight option; that copy and its test
`heatNd_test_posw` are retired (moved out of the repository). `test_copquadvar_faces`
replaces the test. **Measured**: all 24 SDPs of the switch-over set are bit-identical
before and after (§9), so no number in this README changed.

The benchmark follows D. S. Jagt and M. M. Peet, *A State-Space Representation of
Linear Multivariate PDEs and Stability Analysis using SDP*, arXiv:2508.14840v4 (16 Sep
2026). The paper's own numerical example (Ex. 2/21, Sec. 7.2.1, Table 1) is the N = 2
member of the family below. Stock PIETOOLS stops at two spatial variables (paper
Sec. 7.1; `lpiprogram.m`). This benchmark uses the n-variate container classes
(`copvar`/`cdopvar`) on branch `ndopvar` to pose the same LPI in 1, 2 and 3 variables.

Every number below was **measured** on this workstation unless it is marked
*inferred*. Setup: MATLAB R2025b with `maxNumCompThreads = 8`, and MOSEK 11.0.30.
MOSEK does not inherit MATLAB's thread limit: it reported 24 threads in every solve of
this revision, recorded per trace row. The implementer's earlier 3-D traces did not
record it; its logs say 24.

## 0. Library state the results depend on

**Moved and committed (CC, 09/28/2026).** The folder moved from
`sopvar/Testfolder/sdopvar/claude_tests/heatNd` to `PIETOOLS_demos/sopvar_demos` and is now
tracked; `heatNd_path` finds the root two levels up. Checked with `git` on 09/28/2026 (HEAD
`b3334203`): `Sedumi2Mosek.m` (a035aa8e) and `executives/*.m` (e6c78800) are committed, and
`sossolve.m` and `repmat.m` have no working-tree changes. `copquadvar.m`, `lpi_eq_cdopvar.m` and
`lpi_eq_sdopvar.m` are committed in e792cc40 (09/28/2026). Before that commit, `copquadvar`
rejected the face codes ('psatz' should be 0 or 1), so every `'linear'` / `'faces'` Psatz term
errored: `heatNd_lpi` (default preset `'bench'`), `heatNd_bisect`, `heatNd_validate`,
`heatNd_ladder`, `heatNd_poincare` (default gen `'faces'`), `heatNd_poincare_validate` (all SDP
stages but `'3d0'`), and `test_copquadvar_faces`, `test_heatNd_lpi`, `test_heatNd_poincare`
(read from the code). Measured on e792cc40 with a clean path: those three tests pass (97, 74,
44 checks), and `cx_run_1d` (21 stock-vs-container 1-D executive cases) prints identical
verdicts, brackets, shapes and residuals with and without the batched `soseq`. The 09/27
record follows unchanged.

Measured with `git` on 09/27/2026: branch `ndopvar`, HEAD `5234b627`. The heatNd folder
is untracked. These uncommitted working-tree files are on the benchmark's call path:

- Functional changes:
  - `SOSTOOLS400/internal/processing/Sedumi2Mosek.m`: the O(nnz) `bara` construction
    (MMP 09/26/2026). `heatNd_solve` calls it directly. At HEAD, the per-row `cellfun`
    conversion would cost minutes per 3-D solve (*inferred* by review, not run).
  - `sopvar/lpis_sopvar/lpi_eq_cdopvar.m` and `lpi_eq_sdopvar.m`: batched `soseq`
    (MMP 09/26/2026). Every equality of `heatNd_lpi` goes through them.
  - `sopvar/lpis_sopvar/copquadvar.m`: the face codes `options.psatz` = 2d+1 / 2d+2
    (MMP 09/27/2026). Every `'linear'` / `'faces'` Psatz term of `heatNd_lpi` and
    `heatNd_poincare` uses them (§9). With psatz 0 or 1, its outputs are bit-identical to
    the committed version (measured, §9).
- Comment-only changes (0 non-comment lines changed): `SOSTOOLS400/sossolve.m` (the b
  normalization that `heatNd_sdp` replicates is committed) and
  `SOSTOOLS400/multipoly/@polynomial/repmat.m`.

Also uncommitted in the tree, but NOT called by heatNd (checked with grep):
`executives/*.m` (7 files), `CLAUDE.md`, and `PIETOOLS_demos/cuadmm/*` and
`PIETOOLS_demos/lowrank_2d_stability/*`. Commit or record the listed files before
publishing numbers from this benchmark.

## 1. PDE family and exact answers

    u_t = sum_{i=1}^N u_{s_i s_i} + r u   on [0,1]^N,
    s_1: u = 0 at 0 and 1 (DD);   s_2..s_N: u = 0 at 0, u_s = 0 at 1 (DN).

`heatNd_pie` also accepts any per-direction BC from {DD, DN, ND} and any box
`prod [a_i,b_i]`. NN is excluded because d²/ds² is not invertible there.

Separation of variables gives the following (paper App. C.1, extended to N directions):

- λ_i = π²/L_i² for DD and π²/(4L_i²) for DN or ND.
- λ1 = Σ λ_i, which is the lowest eigenvalue of −Δ.
- The PDE is stable iff r < λ1, and the exact decay rate is **k\* = λ1 − r**.
- Eigenfunction: sin(πs_1) Π_{i>1} sin(πs_i/2).
- Gain: ‖u(t)‖ ≤ (1/Πλ_i) e^{−k\* t} ‖D^δ u0‖.

### 1.1 r only shifts the rate: the benchmark parameter is κ = r + k

A = A0 + rT (§2), so the LPI of Cor. 35 at (r, k) is the LPI at (0, κ) with κ = r + k:

    P*A + A*P + 2k P*T = P*A0 + A0*P + 2(r + k) P*T,

and (15a) does not involve r. The SDP depends on κ only. Measured:

- The reviews found bitwise-identical SDPs at (r, k) and (0, r + k) in 1-D, 2-D and 3-D.
- `heatNd_lpi` now builds (15b) from A0 at κ = r + k. It reproduces the earlier dumps
  bitwise in 1-D and 2-D, including an r = 12 dump rebuilt at r = 0 and κ.
- `test_heatNd_lpi` (7) composes the literal Cor. 35 (15b) with A = rT + A0 and compares
  it with the κ build: equal to 0 at r = 12 and to 2.8e−17 relative at r = 3.7. It also
  compares the κ build with the r = 0 instance built from its own base: equal to 0.

So the benchmark is parameterized by (N, preset, d) and κ. The threshold is κ\* = λ1,
and **the r-invariant reach is the gap k\* − k̂ = λ1 − κ̂**. An r ≠ 0 is accepted
everywhere, but it only relabels k = κ − r. Relative reach k̂/k\* depends on r and is not
used. Thm. 34 needs k ≥ 0, so at a given r only κ ≥ r is a rate. For r ≥ λ1 (k\* ≤ 0),
`heatNd_bisect` raises an error; bisect the r = 0 instance instead, since it is the same
SDP.

| N | BC | λ1 = κ\* |
|---|---|---|
| 1 | DD | π² = 9.869604 |
| 2 | DD×DN (paper Ex. 2) | 5π²/4 = 12.337006 |
| 3 | DD×DN×DN | 3π²/2 = 14.804407 |

## 2. PIE (paper Thm. 19)

The PIE state is v = D^δ u = u_{s1s1…sNsN}. It evolves by T v_t = A v, with:

- T = T_N⋯T_1, where T_i = (d²/ds_i²)^{-1} on the lifted 1-D domain.
- A = r T + A0, with A0 = Σ_i Π_{j≠i} T_j, because d²/ds_i² T_i = I.

T_i has the Green kernel G = h + (s−θ)[θ ≤ s] (paper Cor. 4/20), where:

- DD: h = −(s−a)(b−θ)/L.
- DN: h = a − s.
- ND: h = θ − b.

On [0,1] these are the paper's Ex. 3 kernels. Because the T_i act in disjoint variables,
the kernel of T is the product of the 1-D kernels. `heatNd_pie` therefore builds T, A
and A0 directly as 1×1 `copvar` objects on L2[s1..sN], with
`params{g1..gN} = kron(C^1_{g1},…,C^N_{gN})`. Only 2^N of the 3^N cells of T are nonzero.

## 3. The LPI posed (paper Cor. 35) and its encoding

Find P = c'Z_d (indefinite) and R, Q ⪰ 0 such that

    (15a)  P*T − ε² T*T − R = 0
    (15b)  P*A + A*P + 2k P*T + Q = 0.

If these hold, the PDE is exponentially stable with rate k (Thm. 34, which needs ε > 0
and k ≥ 0; `heatNd_lpi` rejects ε ≤ 0 and k < 0). The paper uses ε = 0.1.

`heatNd_lpi(pie,d,k,ep,opts)` returns the **unsolved** program at κ = pie.r + k, with
`eq` expressions only.

**Encoding.**

- (15a) is imposed in the Sec. 7.1 listing's form, `'split'`: P*T − T*P = 0 and
  P*T − ε²T*T − R = 0, each on one cell per adjoint pair. Measured against the
  every-coefficient form `'full'`, it has the same row rank, and the stacked systems have
  the same rank, so the affine set is the same (N = 1, 2 here; N = 3 by review). The row
  count depends on the bases. With the default `'Tspan'` bases split has 20–30% fewer
  rows, but with `'paper'` dp = 0 at d = 1 it has **more** (review: 2-D 277 → 301, 3-D
  3593 → 4295). The verdicts matched at every k tested.
- (15b) is X1 + κX2 + Q = 0 with X1 = P*A0 + (P*A0)* and X2 = P*T + (P*T)*, imposed
  `'symmetric'`.
- The SDP is **affine in κ**: At(κ) = At(κ1) + (κ − κ1)·dA, with b and K fixed. This was
  measured exact to 0 in 1-D and 2-D, and in 3-D by review. `heatNd_bisect` relies on it.

**Bases.**

- **P.** P = c'Z_d uses exactly the paper's Z_d (eqs. 13–14): joint degree ≤ d in (s,θ)
  per direction, giving μ(d)^N decision variables with μ(d) = d²+4d+3. `'tensor'`
  instead uses `lpivar_cdopvar`, which is the stock `lpivar` convention.
- **R and Q.** These use Z_{d'}* M[·] Z_{d'}. The default basis, `'Tspan'`, takes per
  direction the smallest set {θ^a s^c : a ≤ 1, c ≤ 1, a+c ≤ c_j} whose span contains
  T_i: c_j = 2 for DD and 1 for DN/ND, read off T_i's kernels. Every cap is raised by
  `dp`. T is then in the span of R's basis, so T*T is a Gram form in it. P = T itself
  is in the P space only for d ≥ 2 in a DD direction and d ≥ 1 in a DN/ND direction,
  because the DD kernel sθ has joint degree 2. Other choices:
  - `'paper'`: Z_{d'} exactly.
  - `'tensor'`: no joint cap.
  - `'balance'`: the listing's `degbalance` rule.

  Only the basis blocks that can reach the cancelled operator's support are declared
  (`eq_opts_sopvar`). That pruning is argued lossless in that file. Its one measured
  comparison (review, 2-D bench d = 0) certified *less* with pruning off (12.045447)
  than on (12.182793), although pruning off has the larger feasible set. That was the
  old, solver-limited bisection (§5) and does not measure the pruning. Re-measured with
  the current bisection (tight, retry, rtol 1e−5): pruning off gives [12.183790,
  12.184449), with m 1400 and Gram 10×80 and 3 uncertain points inside; pruning on
  gives [12.183790, 12.183978). The certified lo is the same, so this case shows no loss
  from pruning. Losslessness in general is still only argued.

**Psatz terms.** These are **not in Cor. 35 as printed.** Positivity is required only on
the domain, through R = Z*Q0Z + Σ_g Z* g Q_g Z with g ≥ 0 on the box. `'linear'` (the
default) uses the 2N normalised face generators (θ_i−a_i)/L_i and (b_i−θ_i)/L_i. This is
stock `poslpivar_2d` psatz = [3 4 5 6] (CC, 09/23/2026), extended to N-D as `poscopvar`
psatz = 2i+1 (lower face of s_i) and 2i+2 (upper face), one positive variable per face.
`'product'` is stock psatz = 1, and `'none'` is the paper as printed.

These are the presets. The reach is given as κ̂ (= k̂ at r = 0) and the gap λ1 − κ̂,
from §6:

| preset | R, Q basis | Psatz | 1-D gap | 2-D gap | 3-D |
|---|---|---|---|---|---|
| `bench` (default) | Tspan, dp = 0 | linear | 0.1531 (d = 0, 1) | 0.1530–0.1532 (d = 0, 1) | d = 0: gap in (0.104, 0.204] (κ̂ ∈ [14.60, 14.70)) |
| `heavy` | Tspan, dp = 1 | linear | ≤ 5.2e−5 (d = 0, 1) | ≤ 6.5e−5 (d = 0, 1) | no: build priced at 450 GB, skipped |
| `listing` | balance | none | none certified | none certified | – |
| `listingL` | balance | linear | – | 0.160 (d = 0), ≤ 1.6e−4 (d = 1) | – |

**Why Psatz is needed (measured).** Cor. 35 as printed, with no Psatz term, certifies no
positive rate in 1-D or 2-D at the degrees tried (§6). This matches the scout's stock
reproduction, where the verbatim listing reproduces no row of Table 1. In 2-D, from the
review's ablation of bench d = 0:

- The product generator is certified infeasible at 0.5 κ\*.
- Dropping the Neumann face, or the s2 Dirichlet face, gives **no certificate** (MOSEK
  UNKNOWN, rel_b 3.5e−3 to 7e−3, no valid ray) at κ = 0, 0.1 κ\* and 0.5 κ\*. This is not
  a proof that the face is needed.
- Putting the Psatz term on Q only (`psatzon 'Q'`) is **certified infeasible** at the
  same κ (Farkas relative violation ≤ 9e−13, under the earlier relative gate of §5; the
  absolute violation was not recorded). Psatz on R is needed for feasibility, not
  just for accuracy.
- Dropping P's multiplier cells is certified infeasible at 0.5 κ\*.

## 4. Files and how to run

Paths are set up by `heatNd_path` (below); every other entry point calls it itself.
Run each entry point from a clean working directory (for example this folder), with a
command such as:

    matlab -sd "<this folder>" -batch "test_heatNd_pie"

Semantics and the building blocks:

| file | what it does |
|---|---|
| `heatNd_path.m` | Puts this checkout on the path and removes `.claude`/`worktrees` entries. It asserts that every routine the benchmark relies on resolves once, inside the checkout. |
| `heatNd_pie.m` | `pie = heatNd_pie(N,r,bc,dom)` returns T, A and A0 (1×1 `copvar`), the vars, the domain, and `exact` (λ, λ1, k\*, eigenfunction, T/A eigenvalue factors). |
| `heatNd_apply.m` | Applies a fixed `sopvar`/`copvar` to a function by quadrature, from the class definition. This is the tests' evaluator. |
| `heatNd_lpi.m` | `[prog,meta] = heatNd_lpi(pie,d,k,ep,opts)` returns the unsolved Cor. 35 program at κ = pie.r + k (`meta.kappa`). `meta` holds stage times, memory, q, nP and the SDP shape. Passing `opts.base` reuses the κ-independent part (also across r). |
| `heatNd_sdp.m` | `D = heatNd_sdp(prog,file,meta)` returns the pre-solve SeDuMi SDP (At, b, c, K) that `sossolve` would build, with nothing shadowed. It saves that SDP, or a given SDP struct, with metadata (v7.3): κ, the independent row set and the reference verdict with its route (§8). |
| `heatNd_solve.m` | `v = heatNd_solve(prog, D or file, opts)` solves with MOSEK and returns a **certified** verdict (§5). `tight` sets MOSEK tolerances to 1e−10. For a dump file, the stored row set and reference route are the defaults. `route 'lpisolve'` cross-checks through `sossolve`. |
| `heatNd_bisect.m` | `B = heatNd_bisect(pie,d,ep,opts,bopts)` bisects on κ over certified verdicts (§5), using the κ-affine SDP. It can resume (`trace0`, `kaff_file`), and it accepts an `oracle` in place of the solver for testing. |
| `heatNd_ladder.m` | `L = heatNd_ladder(rungs,lopts)` is a build-only sweep over (N, d, preset). It **prices each rung before building** and skips any projected over 24 GB or 20 min, logging what it skipped. Optional dumps (off by default) are named by κ and carry no reference verdict. |
| `heatNd_lindep.m` | `keep = heatNd_lindep(At)` returns a maximal independent set of equality rows, from one sparse QR. It matched the SVD rank in every 1-D and 2-D case tested. MOSEK's presolve removes none of the dependent rows. |
| `heatNd_validate.m` | `heatNd_validate('1d'/'2d'/'3d',vopts)` runs the validation below and writes straddling dumps. `vopts.dir` and `vopts.dumpdir` default under `tempdir`, which is **not durable**. `'3d'` resumes from a saved trace and can take a list of κ (`klist`). |
| `heatNd_poincare_ops.m` | The Poincaré diagnostic's operators A_i, A0 and D_iT, and the SDP-free identity checks I0 (§11.1). |
| `heatNd_poincare.m` | The λ-affine Poincaré LPI family (L2 and energy forms, full or per direction); modes build / test / bisect / objective / sdp (§11.2). Its bisection and tests run through `heatNd_bisect` (solve hook `bopts.oracle`). |
| `heatNd_poincare_validate.m` | The Poincaré validation stages `'I0'`, `'1d'`, `'2d'`, `'2d2'`, `'3d0'`, `'3d'` and `'3dsum'` (§11.3–11.4). |

Tests. Each one asserts its checks and errors if any fails.

| test | checks |
|---|---|
| `test_heatNd_pie.m` | For N = 1, 2, 3, including non-default BCs and a non-unit box, 70 checks against each object's own definition: T v equals the exact solution of D^δ u = v with the BCs; T D^δ φ = φ; the BC traces; A v and A0 v; the eigen identities T φ = Π(−1/λ_i) φ and A φ = (r−λ1) T φ; the paper's Ex. 3 kernel; and STOCK `convert(PDE,'pie')` T and A for N = 1, 2, evaluated from the opvar/opvar2d definitions with no converter. It also has a negative control. |
| `test_copquadvar_faces.m` | The face codes of `copquadvar`/`poscopvar` that the benchmark uses (§9). (i) Value checks for N = 1, 2, 3: psatz 2, 2N+3, 1.5, [3 4] and type `'sym'` with a face code error, and 0, 1, 3..2N+2 are accepted. (ii) ⟨x,Wx⟩ = ∫ g (Zx)'Q(Zx) at degree 0 by independent quadrature, with g built from the definition, for N = 1, 2, 3 on a non-unit box and every code. (iii) A control: the other face of the same direction must not match. (iv) The Galerkin matrix on tensor monomials is symmetric PSD (N = 1, 2). It replaces `heatNd_test_posw`, which tested the retired copy. `test_heatNd_lpi` (6) runs it. |
| `test_heatNd_lpi.m` | 74 checks: basis counts against the paper formulas; the split/full rank and verdicts; affinity in κ; that the dump route agrees with `lpisolve`; that **a certified point is a semantic certificate** (P, R, Q at the returned values satisfy the operator inequalities on test functions and on the lowest mode); (7) the literal Cor. 35 against the κ build, and against the r = 0 instance from its own base; (8) the bisection's bracket logic on synthetic verdicts (an uncertain hole does not cap the search, an uncertain band is never "resolved", the retry is taken, the λ1 cap holds with its exact stop label, a certified +1 above λ1 is rejected and stops the run, the upward doubling stops at `kmax`) and the input guards; (9) the Farkas gate on the absolute violation; (10) a dump file repeats its stored route, and a given `keep` overrides it. |
| `test_heatNd_poincare.m` | 44 checks of the Poincaré diagnostic (§11): I0 for N = 1, 2, 3; three negative controls that must fail in the expected I0 layer; SDP shapes against the pricing formulas and stock; op = grad; the d = 1 face values; soundness at 1.001 × exact; coverage at d = 0; the bisection through `heatNd_bisect`; and semantic certificates: the Galerkin matrix of G at the solution against the analytic form, and its Rayleigh–Ritz minimum over the tensor polynomial space (> 0), with the analytic form at 1.001 × exact as a control that must be negative there. |

## 5. Verdict rule and bisection

**Verdicts (never solver status alone).** The solve is certified feasible (+1) only if
all of the following hold:

- MOSEK returns PRIMAL_AND_DUAL_FEASIBLE / OPTIMAL.
- rel_b = ‖At'x − b‖/‖b‖ ≤ 1e−6 on **all** rows of the dumped SDP.
- Every PSD block has min eig ≥ −1e−8 × (the largest eigenvalue over all blocks).
- x ≠ 0.

It is certified infeasible (−1) if either of these holds:

- MOSEK returns PRIMAL_INFEASIBLE_CER, **and the ray is verified**. Scaled to b'y = −1
  (‖b‖ = 1), the **absolute** violation of At*y ∈ K\* must be ≤ 1e−8; otherwise the
  verdict is 0, "not verified". With w = At*y, a verified ray excludes every x with
  ‖x_f‖ + Σ_j tr(X_j) below `cert_radius` = 1/max(‖w_f‖, the largest negative eigenvalue
  of the W_j). The gate used before divided the violation by the ray's own slack scale.
  The final review showed that a crafted ray with violation 1 on a feasible SDP passes
  that test, so it was replaced. Measured by the review:
  - the certified −1s of the recorded 1-D/2-D traces have absolute violations 1.2e−14 to
    3.8e−10 (relative 2.7e−18 to 6.1e−13);
  - re-solves of every 1-D/2-D dump give ≤ 1.9e−9, i.e. radii 6.6e8 to 4.2e12, against
    tr(X) of 3.5e2 to 6.9e3 at certified +1 points (1.58e4 in 3-D);
  - the 3-D −1s have 1.4e−12 (κ\*) and 6.1e−12 (14.70).

  So the new gate changes no recorded verdict.
- The SDP has a row with no variable and b ≠ 0. This is flagged as `deficient`: a degree
  deficiency.

Anything else is uncertain (0). It is reported and never coerced.

**Accuracy floor (measured).** With MOSEK's default tolerances, rel_b on these SDPs sits
at 1e−7 to 3e−6, right at the 1e−6 gate. The reviews saw verdicts at the same κ flip
between 0 and +1 under 1e−15 changes of the data. With `tight` (MSK_DPAR_INTPNT_CO_TOL_
PFEAS/DFEAS/REL_GAP = 1e−10), rel_b falls to 1e−10 to 1.2e−8 at 8 points near κ̂ in 2-D
(bench, heavy; all rows and independent rows). That costs about 5 more iterations, or
+35% MOSEK time. Example, 2-D heavy d = 1 at κ = 12.3274: defaults on all rows give 0
(rel_b 1.55e−6), while tight gives +1 (8.9e−9). Near κ\*, MOSEK often stalls at the same
iterate under either tolerance and returns UNKNOWN. In the recorded runs, just above λ1
that happened on both row sets.

**The +1 gate can accept infeasible SDPs within about 1e−8 of the threshold** (measured
by the review, tight tolerances, direct builds): 1-D heavy d = 0 is +1 at λ1 + 1e−9
(rel_b 2.5e−8), and 2-D heavy d = 0 is +1 at λ1 + 1e−8 (rel_b 4.1e−7, OPTIMAL). rel_b
alone does not separate these from genuine certificates: the 3-D +1 at 14.60 has 3.3e−7.
The false-acceptance width in κ is about 1e−8 to 1e−7 (*inferred*), far below every
reported bracket. The bisection rejects such points (Cap, below).

**Bisection (on κ, certified verdicts only).**

- lo is the largest κ with a certified +1, and hi the smallest with a certified −1.
- An uncertain point **never bounds the search**. The next trial is the midpoint of the
  widest gap between consecutive points of {lo, uncertain points in (lo, top), top}.
- "Resolved" only when a certified −1 closes the bracket: hi − lo ≤ max(rtol·|lo|,
  1e−6). Uncertain points are reported separately (`unc` in (lo, hi), a missing end read
  as ±∞; `unc_all`), and `stop` says why the search ended.
- **Retry.** An uncertain verdict is re-solved per `retry`, cumulatively: `'tight'`,
  `'loose'` (MOSEK defaults) and `'rows'` (the other row set). Validation uses tight
  tolerances throughout, with retry {rows, loose}. Every attempt is a trace row.
- **Cap.** Trial points stay at or below min(hi, λ1 + tol). Every κ > λ1 is infeasible,
  since a certificate would give a rate above the exact k\* (Thm. 34). In the recorded
  runs MOSEK returned UNKNOWN there (2-D listingL d = 1 spent its 30-attempt budget on 8
  such points).
  - A certified +1 above λ1 is therefore a contradiction, i.e. a numerical acceptance by
    the gate. It is kept out of lo and recorded in `above_kstar`, and the run stops with
    "certified +1 above lambda1: numerical acceptance" (not resolved). Before this fix, a
    review run (1-D heavy d = 0, rtol 1e−9) took a +1 at λ1 + 6.8e−10 as lo, under a
    benign stop label and with a negative gap(2).
  - lo = λ1 exactly stops with "lo = lambda1 (Thm. 34 bound)". gap(2) is clipped at 0,
    and `gap_clipped` flags the clip (only possible with the cap off).
  - The verdicts stay certified; the theorem only decides where to solve and which +1 to
    reject.

  Consequently, when κ̂ is within tol of λ1, "resolved" needs a certified −1 within tol
  above lo, which MOSEK rarely gives. The operative result is then the gap bound
  [0, λ1 − lo], where the 0 comes from the theorem.
- **Upward doubling.** Without a certified −1, κ is doubled up to `kmax` (default 4 λ1),
  then the run stops with "no certified -1 up to kmax". Unbounded, an UNKNOWN band above λ1
  spent the budget at κ up to 1.1e4 (review, oracle).
- Every trace row records κ, verdict, rel_b, psd_relmin, MOSEK time, iterations, η,
  **MOSEK's thread count**, the absolute Farkas violation (the gated quantity), the row
  set, tight, and the attempt number.

The earlier bisection let the smallest uncertain point above lo cap the search, and could
report "resolved" on it. Its reported k̂ values mixed LPI reach with MOSEK accuracy. §6
shows which values moved.

## 6. Measured validation (r = 0, so k = κ; ε = 0.1)

Runs: `heatNd_validate('1d')` and `('2d')` with tight tolerances, retry {rows, loose},
rtol 1e−5 and tmax 300 s (480 s for 2-D heavy) per bisection. 24 MOSEK threads in every
solve. Brackets are [lo, hi): lo is certified +1 and hi certified −1. The gap is
λ1 − κ̂ ∈ [max(0, λ1 − hi), λ1 − lo]. The "earlier" column is the implementer's
default-tolerance bisection at rtol 1e−3, converted to κ = r + k; values marked † were
capped by an uncertain point.

### 6.1 1-D (N = 1, DD, κ\* = 9.869604)

| preset | d | κ bracket (this run) | gap | stop | earlier, r = 0 / r = 9 as κ |
|---|---|---|---|---|---|
| bench | 0 | [9.716459, 9.716534) | 0.1531 | resolved | 9.711537 [9.717319] / 9.716065 [9.716574] |
| bench | 1 | [9.716459, 9.716534) | 0.1531 | resolved | same as d = 0 |
| heavy | 0, 1 | lo 9.869552 (hi 10.856565; unc 9.869628 > λ1) | ≤ 5.2e−5 | cap: uncertain above λ1 | 9.867677 [9.873460] / 9.869435 [9.869944] |
| product (listing + product Psatz) | 0 | [7.999964, 8.000040) | 1.870 | resolved | 7.999777 [8.005560] / none (κ = 9, k = 0, certified −1) |
| product | 1 | [9.857579, 9.857880) | 0.0117–0.0120 | budget (3 uncertain inside) | 9.856111 [9.861894] / 9.857206† [9.858734] |
| productT (product, tensor P = stock analogue) | 0 | [7.999964, 8.000040) | 1.870 | resolved | 7.999777 [8.005560] / none (κ = 9 certified −1) |
| productT | 1 | [9.869552, 9.869628) | ≤ 5.2e−5 | resolved | 9.867677 [9.873460] / 9.869435 [9.869944] |
| listing | 0, 1 | none; certified −1 at 0.019277; κ = 0 uncertain | – | no feasible κ | the same |

A first run with the same code but without the λ1 cap certified heavy's −1 at 9.869665
(bracket [9.869574, 9.869665), resolved). The two runs' certified points are mutually
consistent.

These are consistent with the stock brackets from the scout ref, on this tree:

- The stock verbatim listing gives none at d = 0 (k = 0 uncertain) and [0, 0.0193] at
  d = 1. The container gets the **same** infeasible point, 0.019277.
- Stock listing with psatz = 1 gives [7.99995, 8.00010] at d = 0 and
  [9.869526, 9.879604] at d = 1 (at r = 9 as κ: none at d = 0, [9.869526, 9.879604] at
  d = 1, the same SDPs).

No recorded run certified any κ above κ\*. That describes these grids, not the procedure:
the gate can accept an infeasible SDP within about 1e−8 of the threshold (§5), and with the
cap such a point now stops the run instead of becoming lo.

### 6.2 2-D (N = 2, DD×DN, the paper's example, κ\* = 12.337006)

| preset | d | κ bracket (this run) | gap | stop | earlier, r = 0 / r = 12 as κ |
|---|---|---|---|---|---|
| bench | 0 | [12.183790, 12.183978) | 0.1530–0.1532 | 1 uncertain inside | 12.182793 [12.190022] / 12.183510 [12.184300] |
| bench | 1 | [12.183790, 12.183978) | 0.1530–0.1532 | 1 uncertain inside | the same as d = 0 |
| listingL | 0 | [12.176354, 12.177107) | 0.1599–0.1607 | budget (4 uncertain inside) | **12.030990†** / 12.172452† |
| listingL | 1 | lo 12.336846 (hi 13.570706; unc 12.336941, 12.337035) | ≤ 1.6e−4 | cap | **12.327367†** / 12.331213† |
| listingLT (listingL, tensor P) | 0 | as listingL d = 0 | 0.1599–0.1607 | budget | 12.030990† / 12.172452† |
| listingLT | 1 | as listingL d = 1 | ≤ 1.6e−4 | budget | 12.327367† / 12.281452† |
| heavy | 0 | lo 12.336941 (hi 13.570706; unc 12.337035 > λ1) | ≤ 6.5e−5 | cap: uncertain above λ1 | **12.334596†** / **12.318576†** |
| heavy | 1 | lo 12.336941 (hi 13.570706; unc 12.337035 > λ1) | ≤ 6.5e−5 | cap: uncertain above λ1 | **12.320139†** / **12.315416†** |
| listing | 0 | none; certified −1 at 0.048191; κ = 0 uncertain | – | no feasible κ | none; k = 0 uncertain |
| listing | 1 | none; certified −1 at 0.385531 | – | no feasible κ | none |

A first run of the same code without the cap gave bench [12.183809, 12.184035) and
listingL d = 0 [12.176468, 12.177033). Its certified points are consistent with those
above, so the union bracket for bench is [12.183809, 12.183978).

**Which earlier values were solver-limited** (the earlier bisection stopped at an
uncertain point, and the value moved by more than its tolerance):

- listingL d = 0: 12.030990 → 12.176354 (+0.145).
- listingL d = 1: 12.327367 → 12.336846.
- heavy d = 0: 12.334596 → 12.336941, and heavy d = 1: 12.320139 → 12.336941. The
  earlier r = 12 runs (as κ, 12.318576 and 12.315416) were lower still. The reviews
  had already certified 12.335096 and 12.3346 on direct builds. The earlier "heavy
  d = 1 below d = 0", though the d = 1 P space contains d = 0's, was solver noise.
  At d = 1, the point 12.336941 was 0 with tight tolerances on all rows (rel_b 1.0e−5)
  and +1 on the independent rows (3.6e−8): the row retry decided it.
- The earlier r = 0 and r = 12 runs are the same SDPs, and they disagreed by up to
  0.14 (listingL d = 0). That spread is solver noise, not an r effect.
- bench was not solver-limited. It moved only within the earlier rtol-1e−3 bracket.

d = 0 and d = 1 give the **same** result for bench and heavy, in 1-D and 2-D. Both use
the fixed `Tspan` R, Q basis. They do not for listingL (gap 0.16 against ≤ 1.6e−4) or
for 1-D product (1.87 against 0.012). Those use the `balance` basis, which is sized from
the operator and grows with d. The link to the basis is *inferred* from these five
cases.

Comparison with the paper's Table 1 (MOSEK, i7-5960X), converted to κ = r + k:

| | d = 0 | d = 1 | d = 2 |
|---|---|---|---|
| r = 0 | 12.336 | 12.336 | 12.337 |
| r = 12 | 12.33578 | 12.33666 | 12.33700 |

The paper's own r columns show the κ invariance: k̂ + r = 12.336, 12.3362, 12.3340,
12.33578 and 12.33613 across its r values (review). heavy certifies κ = 12.336941 at
d = 0 and d = 1. That is above every Table 1 entry for d = 0 and d = 1, and below the
d = 2 entries (12.337 and 12.33700) by at most 6e−5, which is within heavy's gap
bound. listingL d = 1 certifies 12.336846, above Table 1's d = 1 entries. The
earlier README's "heavy reaches Table 1 to 5% at r = 12" compared
solver-limited values as relative reach, and is withdrawn. Neither the stock
reproduction nor this container reproduces Table 1's d = 0 and d = 1 rows **with the
listing as printed** (no Psatz term). Table 1 was computed with an earlier version of
PIETOOLS, whose settings need not carry over to this tree (maintainer, 09/27/2026), so
this is a version difference, not a defect of the benchmark (§10).

Comparison with stock on this tree (scout ref, as κ):

- The verbatim listing is infeasible at k = 0 for both r.
- The listing with linear psatz gives, at r = 0: [10.124, 11.062] (d = 0) and
  [12.084, 12.297] (d = 1). At r = 12 it is infeasible at 0 (d = 0) and gives
  [12.13345, 12.29701] (d = 1).

The stock r = 0 and r = 12 brackets are mutually consistent as κ: κ = 12 lies above
the d = 0 bracket, and the d = 1 brackets overlap. They came from a default-tolerance
bisection, so they may be solver-limited as the container's earlier values were
(*inferred*; not re-measured with tight tolerances).

No recorded run certified any κ above κ\* (a statement about these grids; see §5 and §6.1).

### 6.3 3-D (N = 3, DD×DN×DN, κ\* = 14.804407), bench, d = 0

The SDP has m = 31770 rows, 2.7M decision variables and 7×288 + 7×552 Gram blocks (§7).
`heatNd_lindep` finds rank 22226. That costs one sparse QR of 363 s, done once and
stored with the κ-affine SDP. MOSEK then sees only the independent rows. **Every verdict
is computed on all 31770 rows.**

**Merged κ trace.** The implementer's r = 0 and r = 12 bisections are one SDP family,
merged here with κ = r + k. All were solved on the independent rows with MOSEK defaults,
except the first κ = 11.84 solve (all rows). Thread count was not recorded per solve then
(24 per its log):

| κ | κ/κ\* | source | verdict | residuals | MOSEK s |
|---|---|---|---|---|---|
| 7.402203 | 0.500 | r = 0 | +1 | rel_b 6.8e−7 | 440 |
| 9.622865 | 0.650 | r = 0 | +1 | rel_b 5.1e−7 | 472 |
| 11.843525 | 0.800 | r = 0 | 0 (twice) | rel_b 1.4e−6 all rows / 1.6e−6 independent rows | 680 / 404 |
| 12.841322 | 0.867 | r = 12 | +1 | rel_b 4.9e−7 | 455 |
| 13.261983 | 0.896 | r = 12 | +1 | rel_b 2.5e−7 | 432 |
| 13.323966 | 0.900 | r = 0 (probe, not in trace) | 0 | MOSEK max violation 3e−6, all rows | 499 |
| 13.682644 | 0.924 | r = 12 | +1 | rel_b 7.1e−7 | 437 |
| 14.064186 | 0.950 | r = 0 | +1 | rel_b 6.6e−7 | 586 |
| 14.243525 | 0.962 | r = 12 | +1 | rel_b 3.0e−7 | 421 |
| 14.434296 | 0.975 | r = 0 | +1 | rel_b 6.9e−7 | 498 |
| 14.523966 | 0.981 | r = 12 | +1 | rel_b 2.7e−7 | 532 |
| 14.804407 | 1.000 | r = 0 and r = 12 | −1 (both) | ray violation 1.4e−12 | 371 / 422 |
| **14.600000** | 0.986 | new, tight | **+1** | rel_b 3.3e−7, η 1.1e−11, psd_relmin +3.3e−14; 25 it, 24 threads, peak 11.45 GB | 697 |
| **14.700000** | 0.993 | new, tight | **−1** | ray violation 6.1e−12 (rel 7.4e−17); 18 it, 24 threads, peak 11.15 GB | 615 |

The uncertain points at κ = 11.84 and 13.32 lie between certified points. With
monotonicity they carry no information (solver accuracy).

**Bracket.** From the stored traces alone, κ̂ ∈ [14.523966, 14.804407), i.e. gap ≤ 0.2804.
The earlier README reported [14.434, 14.804) from the r = 0 trace alone. With the two
new solves (`fix_3d_kap`: the κ-affine r = 0 file, independent rows, tight tolerances,
verdicts on all rows), **κ̂ ∈ [14.60, 14.70)**, that is, 0.986–0.993 κ\* and gap
k\* − k̂ ∈ (0.104, 0.204]. Not resolved to rtol: each further solve costs 10–12 min.

**Gap prediction.** The 1-D and 2-D bench gaps are 0.1531 and 0.1530–0.1532. If the 3-D
gap were the same, κ̂ ≈ 14.651 (*inferred*, from 1-D and 2-D only). The two new solves
test this. The prediction lies inside the measured [14.60, 14.70), so it is
**consistent, to ±0.05**. It would have been refuted by +1 at 14.70 (gap < 0.104) or −1
at 14.60 (gap > 0.204). Nothing here says the gap should be N-independent; the
agreement of 1-D and 2-D is itself only an observation.

**3-D, d = 1 (one earlier solve).** m = 34060, rank 22686 (QR 476 s). At κ = 14.434296
MOSEK (defaults, independent rows) returned OPTIMAL in 612 s with rel_b 5.3e−6, so the
point is **uncertain**. The process's lifetime peak, including the QR, was 22.3 GB. It
was not re-solved with tight tolerances.

## 7. Rung table (measured, `heatNd_ladder`, build only)

These are SDP sizes, independent of r (and so of κ). The earlier table's "r = 12" 3-D
row was the same SDP as r = 0 and is dropped. The columns are:

- q: `numel(decvartable)`, the decision variables.
- nP: the decision variables of P.
- m: the SDP rows.
- K.f: the free variables.
- K.s: the PSD blocks.
- nnz: nnz(At).
- build: the base plus two finalizations.
- MATLAB MB: memory in use after the build.

| N | d | preset | q | nP | m | K.f | K.s | nnz | build (s) | MATLAB MB |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0 | bench | 495 | 3 | 51 | 3 | 8×3, 10×3 | 674 | 1.3 | 2197 |
| 1 | 1 | bench | 500 | 8 | 55 | 8 | 8×3, 10×3 | 741 | 0.5 | 2229 |
| 1 | 2 | bench | 507 | 15 | 62 | 15 | 8×3, 10×3 | 840 | 0.3 | 2237 |
| 1 | 0 | heavy | 1854 | 3 | 87 | 3 | 16×3, 19×3 | 2528 | 0.5 | 2215 |
| 2 | 0 | bench | 40409 | 9 | 1242 | 9 | 48×5, 76×5 | 119495 | 1.9 | 2304 |
| 2 | 1 | bench | 40464 | 64 | 1362 | 64 | 48×5, 76×5 | 122472 | 1.5 | 2306 |
| 2 | 2 | bench | 40625 | 225 | 1714 | 225 | 48×5, 76×5 | 131809 | 1.5 | 2311 |
| 2 | 0 | heavy | 565209 | 9 | 4282 | 9 | 192×5, 276×5 | 1805487 | 5.1 | 2548 |
| 3 | 0 | bench | 2713563 | 27 | 31770 | 27 | 288×7, 552×7 | 17833204 | 69.2 (review: 74.6) | 4006 |
| 3 | 1 | bench | 2714048 | 512 | 34060 | 512 | 288×7, 552×7 | 17935345 | 74.6 | 4062 |
| 3 | 2 | bench | 2716911 | 3375 | 47847 | 3375 | 288×7, 552×7 | 18591090 | 79.1 | 4129 |
| 3 | 0 | heavy | – | – | *~211k* | – | – | *~1.3e9* | – | – |

- The 3-D count was predicted from 1-D and 2-D by the Kronecker model X(3) ≈ X(2)²/X(1)
  before it was built. The prediction for m was 30.2k; the build measured 31770.
- 3-D heavy was skipped: projected 450 GB and 12900 s (*inferred*, Kronecker model).
- Re-measured with the κ build (from A0): 3-D bench d = 0 took 92.1 s for the base plus
  one finalization (R 12.6 s, Q 58.8 s, (15a) 6.3 s, (15b) 12.5 s). MATLAB had
  3757 MB in use and the process peak was 4.3 GB. At, b and K equal the implementer's
  κ-affine file bitwise (max difference 0). So the stored 3-D files are the current
  code's SDPs.

**Solve cost per route** (MOSEK 11, 24 threads). Process peaks are private bytes of the
MATLAB process, which includes MOSEK:

| rung | route | MOSEK time per solve | iterations | process peak |
|---|---|---|---|---|
| 1-D, every preset | all rows | < 0.05 s | 8–40 | – |
| 2-D bench | all rows, tight | 0.4–1.1 s | 19–35 | – |
| 2-D listingL / listingLT | all rows, tight | 1.5–6 s / 8–15 s | 15–30 | 4.2 GB (whole run 'a') |
| 2-D heavy | all rows, defaults / tight | 11–14 s / 16.5–19.6 s | 14–16 / 20–22 | 4.75 GB (probe run) |
| 3-D bench d = 0, m 31770 | all rows, defaults | 499–921 s | 16–25 | **20.1–20.9 GB** |
| 3-D bench d = 0, 22226 rows | independent rows, defaults | 371–666.6 s | 15–25 | **10.7–11.45 GB** (review; the maxima are the §8 dump re-solve) |
| 3-D bench d = 0, 22226 rows | independent rows, tight | 615–697 s | 18–25 | 11.07–11.45 GB (11.07: the final review's re-solve of the 14.60 dump, 617.8 s) |
| 3-D bench d = 1, 22686 rows | independent rows, defaults | 612 s | – | 22.3 GB lifetime, *including* the 476-s QR in the same process |

The independent-row route halves the solve's peak memory and saves 20–40% MOSEK time. The
QR costs 363 s (d = 0) or 476 s (d = 1), once per family. The ladder's solve price
(`sprice`) is calibrated per route on the maxima of these 3-D d = 0 solves: all rows
921 s / 21 GB at m = 31770; independent rows 697 s / 11.45 GB at 22226 rows (it used
586 s / 11.1 GB before the final review, 15–20% below the later keep-route solves). The
m² scaling to other m, and the independent fraction 0.70, are *inferred*. Under that
model, 3-D d = 2 (m 47847) prices at about 2090 s / 48 GB on all rows and 1580 s / 26 GB
on independent rows (*inferred*; the earlier "34 GB" was uncalibrated). It was not solved.

## 8. SDP dumps

**Location.** Dumps go to a directory the caller passes (`vopts.dumpdir`,
`lopts.dumpdir`). The defaults are under `tempdir`, which is **not durable**. The dumps
made so far are in this session's temporary scratch area (ephemeral; copy them somewhere
stable to keep them):

- `…/scratchpad/heatNd/fix/dumps/`: the current dumps, named
  `heatNd_N<N>_d<d>_<preset>_kap<κ>.mat`.
- `…/scratchpad/heatNd/impl/dumps/`: the implementer's 28 dumps, named
  `…_r<r>_k<0.90|1.10>kstar.mat`. They are **superseded**. They carry no row set and no
  route, and their names are misleading: r = 12 files duplicate r = 0 SDPs at κ = r + k,
  the 1-D r = 12 files hold absolute k = 0.9/1.1 (k\* < 0 there), and several pairs do
  not straddle κ̂ (both 2-D bench r = 12 files are −1).

**Format.** v7.3, SeDuMi form: min c'x s.t. At'x = b, x ∈ K.

- `At` (nx × m sparse), `b` (normalized, ‖b‖ = 1), `c` (0), `K` (`f`, `l`, `q`, `s`) and
  `bscl`.
- `info`: a description, m, nx, nnz, `deficient`, `kappa`, the HEATND_LPI meta, and for
  the current dumps:
  - `info.keep`, the independent row set (`heatNd_lindep`);
  - `info.ref`, the **reference verdict and the route that produced it**: st, rows
    ('keep' or 'all'), tight, rel_b or the Farkas violation, MOSEK time, iterations,
    MOSEK threads, and the source.
- `v = heatNd_solve('<file>.mat')` repeats the reference route by default. Pass
  `struct('rows','all')` or `struct('tight',…)` to change it.

**Straddling pairs.** `heatNd_validate` writes, after each bisection, the SDP at lo
(certified +1) and at hi (certified −1), taken from the same κ-affine family the verdicts
were solved on. So each pair straddles the measured κ̂ by construction. That gives 1-D
and 2-D pairs for every preset and d of §6.1–6.2, and single −1 files where nothing was
feasible.

**3-D pairs** (bench d = 0, 53 MB each, with `info.keep` = the 22226 independent rows).
Each file is built from exactly the data its reference verdict was solved on, since
verdicts at the gate flip under 1e−15 changes:

- The tight pair straddling κ̂, written by `heatNd_validate('3d')` from the r = 0
  κ-affine file, independent rows, tight tolerances:
  - `heatNd_N3_d0_bench_kap14.600000.mat`: +1.
  - `heatNd_N3_d0_bench_kap14.700000.mat`: −1.
- The earlier pair (`fix_3d_dump` in the scratch area), MOSEK defaults, independent
  rows:
  - `heatNd_N3_d0_bench_kap14.523966.mat`: +1, from the r = 12 κ-affine file at
    k = 0.9 k\*(12), the implementer's expression.
  - `heatNd_N3_d0_bench_kap14.804407.mat`: −1, from the r = 0 file at k\*.

**Dump route re-measured.** `heatNd_solve('heatNd_N3_d0_bench_kap14.523966.mat')` with no
options took the stored route (independent rows 22226, MOSEK defaults). It returned
**+1**: rel_b 2.686e−7 (reference 2.7e−7), η 1.7e−11, psd_relmin +2.0e−14, 23
iterations. MOSEK took 666.6 s (the reference 531.5 s) on 24 threads, with a process
peak of 11.45 GB. So the documented route reproduces the reference verdict.

The review's re-solve of the old r = 12, 0.9 k\* dump **on all rows** (the only route
that file allowed) gave 0: rel_b 1.155e−6, 921 s, 20.9 GB. The +1 needs the independent
rows, which are now stored in the dump.

## 9. Library notes

The benchmark modifies no library file. The one library change it depends on,
`copquadvar.m`'s face codes, was made by the maintainer's decision (§10) and is
uncommitted (§0).

- **Face codes.** `copquadvar`/`poscopvar` `options.psatz` = 2d+1 weights the Gram form
  by (θ_d−a_d)/L_d, and 2d+2 by (b_d−θ_d)/L_d. Here d indexes the sorted registry, which is
  s1..sN here (`heatNd_pie` allows N ≤ 9). In 2-D with sorted names these are
  `poslpivar_2d` psatz 3–6. The benchmark declares one positive variable per face.
  Until 09/27/2026 it used `heatNd_posw.m`, a copy of `copquadvar` with a `gfun` weight
  option. That copy and `heatNd_test_posw.m` are retired (moved out of the repository
  to the session scratchpad, not deleted). All of the following were **measured**:
  - **Same objects.** `poscopvar(psatz = code)` equals `heatNd_posw(gfun = the benchmark's
    face weight)` byte for byte (prog, Pop, Qcell, basis_list). This holds for
    N = 1, 2, 3, every code, unit and non-unit boxes, degrees 0 and 1 (N = 3: 0 and
    joint 1), and with include/sep.
  - **Same SDPs.** All 24 pre-solve SDPs of the switch-over set are bit-identical
    before and after the switch. The heat LPI set is 13 configs: bench d = 0, 1,
    heavy, listingL, a faces subset, psatzon `'Q'` and non-unit boxes, for N = 1, 2,
    plus the 3-D bench at d = 0 (m 31770, 17.8M nnz). The Poincaré set is 11 families:
    L2, energy, sel, d = 2, a faces subset, the `'op'` target and a non-unit box, plus
    3-D L2 at d = 0 and d = 1 (m 28672, 19.0M nnz). Build times are unchanged (3-D
    heat 60.3 / 60.4 s, 3-D Poincaré 45.1 / 44.6 s).
  - **Defaults unchanged.** With psatz 0/1, `copquadvar`'s outputs are bit-identical to
    the committed version on every call made by `test_poscopvar`,
    `test_poscopvar_vs_poslpivar(_2d)` and `test_poscopvar_stability`, and on the
    container executives' builds.
  - **Same cone as stock in 2-D.** The test is the 2-D L2 Poincaré LPI at d = 1, DD×DN.
    It compares stock `poslpivar_2d` psatz 3–6 (with the container-matched basis of the
    scout harness) against the container codes. This was done for all four faces, for
    faces 4+5, and for faces 3, 5 and 6 alone. In every case m, K.s, nnz(At) and nx are
    equal. The permutation-invariant fingerprints of At and b also agree to 1e−16
    relative: sorted values, per-row and per-column nnz profiles, and row norms. Faces 5
    and 6 differ in nnz (24736 vs 34688), so a lower/upper swap would show. The
    fixed-λ verdicts agree at 1.9–2.1 (full) and 0.9–1.1 (per direction). A 30-attempt
    bisection gives the same certified lo in every case: full 1.9999443, per direction
    0.9999721. The hi ends differ only through uncertain points in the budget-limited
    band above the d = 1 ceiling. Equal fingerprints are strong evidence of the same SDP
    up to permutation, not a proof.
- `lpiprogram` refuses N > 2. `heatNd_lpi` builds what it would return for N = 3, as
  `make_prog` does in `test_lpivar_cdopvar`.
- `lpivar_cdopvar` has no joint-degree cap. `heatNd_lpi` declares the paper's Z_d itself
  (`lpivar_Zd`, with the same layout).
- MOSEK's presolve finds no dependent rows, although the equality systems are
  rank-deficient. Measured in 2-D: rank 990 of 1242 rows (split) and 990 of 1712 (full).
- The 2N face terms in `heatNd_lpi` are summed one `plus` at a time. That is 6.0 s of the
  56.5-s Q stage in 3-D (review); left as is.

## 10. Open questions for the maintainer (not decided here)

- ~~How were the paper's Table 1 d = 0 and d = 1 rows produced?~~ Answered 09/27/2026:
  with an earlier version of PIETOOLS; its settings need not transfer to this tree.
  Neither stock nor the container reproduces those rows here without Psatz.
- ~~Should `copquadvar` get a weight option (retiring `heatNd_posw`)?~~ Decided
  09/27/2026: yes, as the face codes psatz = 2d+1 / 2d+2 (§9). `heatNd_posw` is retired.

## 11. Poincaré diagnostic of the heat LPI (`heatNd_poincare*`)

Added MMP, 09/27/2026. This is a separate, much smaller LPI on the heat LPI's own T and A.
It diagnoses the heat LPI without running it. It follows the scout reference study (the
math, the stock 1-D/2-D numbers, the 3-D pricing and the decision table); every number of
that study used here was re-measured on the container. Run it from this folder like the
rest of the benchmark:

    matlab -sd "<this folder>" -batch "test_heatNd_poincare"
    matlab -sd "<this folder>" -batch "heatNd_poincare_validate('I0')"

The files are listed in §4. `heatNd_poincare` uses `poscopvar` face codes for the face generators
(§9) and `heatNd_sdp`, `heatNd_solve` and `heatNd_bisect` for the verdicts, the certified
bracket and the trace.

### 11.1 Operators and I0 (the first diagnostic; no SDP)

With u = Tv (v = D^δ u):

- A_i = ∂²_{s_i} ∘ T = Π_{j≠i} T_j (identity in s_i); A0 = Σ_i A_i, and A = rT + A0.
- D_iT = ∂_{s_i} ∘ T = (dT_i/ds_i) ⊗ Π_{j≠i} T_j. The 1-D factor has the kernel dG/ds of the
  Green function G on [a,b]:
  - DD: (θ−a)/L below (θ ≤ s), (θ−b)/L above; on [0,1], θ and θ−1;
  - DN: 0 below, −1 above;
  - ND: 1 below, 0 above.
- Integration by parts, with the boundary term zero for DD, DN and ND (one factor vanishes
  at each end): T*A_i = A_i*T = −(D_iT)*(D_iT) ⪯ 0, and T*A0 = −Σ_i (D_iT)*(D_iT).

**I0** (`O = heatNd_poincare_ops(pie,true)`; `O.I0ok`) checks each object against its own
definition:

- (k) The stored kernels of T, A_i and D_iT, in every cell, at random (s, θ), against the
  textbook Green products (G, dG/ds, or the identity in s_i). Tolerance 1e−12.
- (c) Coefficient identities through the class algebra, relative, tolerance 1e−12:
  - A − rT − A0 and A0 − pie.A0;
  - T*A_i − A_i*T and T*A_i + (D_iT)*(D_iT);
  - T*A0 + Σ_i (D_iT)*(D_iT).
- (q) Semantics by quadrature (`heatNd_apply`) against exact polynomial solutions of
  D^δ u = v, tolerance 1e−10:
  - D_iT v = ∂_i u, and A_i v;
  - ⟨v,T*A_i v⟩ = −‖∂_i u‖² and ⟨v,(D_iT)*(D_iT) v⟩ = ‖∂_i u‖²;
  - ⟨v,T*T v⟩ = ‖u‖² and ⟨v,A_i*A_i v⟩ = ‖A_i v‖².

MEASURED: I0 passes at N = 1, 2, 3 for DD, DN and ND, on unit and non-unit boxes (up to 32
checks). The largest error is 3.3e−16 for (k), 2.7e−16 for (c) and 2.3e−15 for (q). It
takes 0.1 s at N = 1 and 5.2 s at N = 3.

A failure implicates the T/A construction both LPIs share: the BC-to-direction map, the
kernel coefficients, the sorted-variable cell indexing, or ZL/ZR (D_iT has numel(ZL) ≠
numel(ZR) in a DD direction, so a ZL/ZR swap is visible). Do not run either LPI until it
passes.

The negative controls in `test_heatNd_poincare` fail as designed:

- Swapped BC labels fail (k) and (c).
- Swapped lower and upper cells of D_1T fail (k), (c) and (q).
- A sign flip of D_2T fails (k) and (q) but **passes (c)**: the identity is quadratic in
  D_iT. That is why I0 has three layers.

### 11.2 The LPIs

For the directions i in `sel` (`'full'` or one index):

- **L2 form**: Q(λ) = Σ_i (D_iT)*(D_iT) − λ T*T ⪰ 0, i.e. ‖∂_sel u‖² ≥ λ‖u‖².
- **Energy form**: Q(λ) = A_sel*A_sel − λ Σ_i (D_iT)*(D_iT) ⪰ 0, i.e.
  ‖Σ_i ∂²_i u‖² ≥ λ‖∂_sel u‖². It is better conditioned: identity plus compact.
- `target 'op'` writes −T*A_sel as −(T*A_sel + A_sel*T)/2. It gives the same SDP as the
  default `'grad'`: At is equal and b agrees to 3e−17 (test (4)).

**Exact constants** (derived in the reference study; the same for both forms):

- Per direction the constant is μ_i, the lowest eigenvalue of −d²/ds_i² under bc{i}: π²/L²
  for DD and π²/(4L²) for DN/ND.
- The full constant is λ1 = Σ_i μ_i.
- Per direction the operators factor as (−T_i − λT_i²) ⊗ T_î² and (I + λT_i) ⊗ T_î². Each
  is ⪰ 0 iff λ ≤ μ_i.
- On [0,1]: μ_DD = π² = 9.869604 and μ_DN = π²/4 = 2.467401.

**Link to the heat LPI** (Cor. 35 at κ = r + k, §1.1):

- With **P = T**, (15a) is (1−ε²)T*T ⪰ 0 and (15b) is 2[T*A0 + κT*T] ⪯ 0. So Cor. 35 with
  P = T certifies **exactly κ = r + k ≤ λ**, where λ is the full L2-form constant: at P = T
  the heat LPI *is* the full L2 Poincaré LPI at λ = κ.
- With **P = −A0** (V = ‖∇u‖²), (15b) is −2[A0*A0 + κT*A0] ⪯ 0, which is exactly the full
  energy form at λ = κ. (15a) is then the L2 form at ε² ≤ λ1, which holds.
- Both certificates are sharp in exact arithmetic. At a finite basis, a heat LPI whose P
  family contains one of them, and whose R and Q bases contain the Poincaré basis and
  generators, certifies at least the Poincaré λ̂ at that basis.
- The heat bench P basis Z_0 contains neither T nor −A0 for N ≥ 2, because T_DD needs
  d ≥ 2. So that inequality is **not** guaranteed for bench; §11.5 says what is measured
  instead.

**Positive variable.** G = Z*MZ + Σ_g Z* g M_g Z with M, M_g ⪰ 0, imposed as G − Q(λ) = 0
by `lpi_eq_cdopvar 'symmetric'`. Both sides are self-adjoint by construction.

- Z is `poscopvar` at tensor degree d (a, c ≤ d per direction, no joint cap).
- L2 form: all-integral multi-indices only ({2,3}^N). This is **lossless in exact
  arithmetic**: the target has no multiplier cell, and a PSD block whose diagonal reaches
  only zero cells is zero. It is not a numerical equivalence (measured by the final
  review). With the delta indices included, the delta blocks have a diagonal of ~4e−7
  rather than 0, and PSD then admits delta × integral cross blocks of size ~√tolerance.
  - They inflate certified values by ~0.1%: 2-D DD×DN L2 full d = 1 faces is +1 at 2.002
    with any delta index, against the pruned λ̂ 2.0000 (the pruned family is 0 at 2.002
    and certified −1 at 2.05). The certificate is still semantically valid, and zeroing
    the delta rows raises rel_b 8537×.
  - They bring the gate to the edge of certifying above μ: 1-D DD L2 faces d = 3 with all
    indices has rel_b 4.9e−7 (< 1e−6) at 1.00001 μ, and is 0 only because MOSEK's status
    was not OPTIMAL. The pruned families had rel_b ≥ 1.2e−5 above μ.

  So the `eq_opts_sopvar` rule must stay on for verdict solves.
- Energy form: all-integral plus exactly one delta in a direction of `sel`. A_i*A_i = I ⊗
  T_î² has the identity in s_i, so all-integral alone would leave its multiplier cell
  uncovered. (This is the reference's 896 = 512 + 3·2·64 count.)
- Both are exactly the `eq_opts_sopvar` support rule on the target, checked on every build
  (`R.include_eqopts`, true in every build).
- Generators:
  - `'none'`;
  - `'product'`: `poscopvar` psatz = 1, Π(θ_i−a_i)(b_i−θ_i);
  - `'faces'`: the 2N linear face weights, `poscopvar` psatz = 2i+1 / 2i+2, i.e. the heat `'linear'`
    Psatz. `spec.faces` restricts them.

  Each generator has its own M_g ⪰ 0 at the full degree d.

**λ enters b only.**

- Each family is built at three λ (0.5, 0.75 and 0.9 of exact) on one positive variable.
- The build asserts that At, K and m are identical and that b is affine. The third-point
  error is ≤ 1e−12; measured 0 to 1.2e−16.
- Every later λ is b(λ) on the same At, recomputing only the `deficient` test.
- At does not depend on λ, so a `heatNd_lindep` row set would be exact for every λ.
  However, `heatNd_lindep`'s sparse QR of the 3-D faces At (1.8M × 28672, nnz 19M)
  exceeded 70 GB within 7 min (MEASURED; killed), so 3-D solves use all rows.
  The systems are rank-deficient: 2-D d = 1 has rank 544/768 (faces) and 356/512 (none);
  3-D d = 1 none has 8896/16384.

**Modes** (`popts.mode`):

- `'test'` (given λ) and `'bisect'` run through `heatNd_bisect` with the family as its
  solve hook (`bopts.oracle`, κ = λ). That gives the certified bracket, retry, cap and
  trace, with MOSEK's thread count per attempt.
- `'objective'`: λ is a new first free variable, with At'x − db·λ = b(0) and max λ. Then
  one fixed-λ re-certification at (1 − 1e−3)λ̂.
- `'sdp'`: the SDPs at given λ.

Near the threshold MOSEK returns UNKNOWN in a band about 1e−4 relative wide, and up to
more than 1e−3 for the 2-D full energy form. There, fixed-λ solves at 0.9999 and 0.99999
× λ̂ = 12.0412727 are 0 under tight and default tolerances, and the first certificate is at
0.999 λ̂ (final review); the certified value is the bisection lo 12.035914. This was
measured in 1-D and 2-D, with tight or default tolerances and on both row sets. A bisection
spends its budget there ("uncertain points fill the bracket"). The objective gets through
the band in one solve, and its optimum is often certified +1 itself; the re-certification
is then a second certificate.

### 11.3 Measured: 1-D and 2-D (`heatNd_poincare_validate('1d'/'2d'/'2d2')`)

All runs used MOSEK 11.0.30 with 24 threads in every solve (MATLAB maxNumCompThreads = 8),
and tight tolerances (1e−10) for every verdict.

Column key:

- λ̂: the objective-mode optimum.
- "cert / re-cert": the verdict at λ̂ itself, and the fixed-λ verdict at 0.999 λ̂.
- "bracket": the certified [lo, hi) of a 40-attempt bisection (retry `'loose'`).
- "stock": the scout reference, bisected to 1e−4 in λ/exact.

Every case was probed at 1.001 × exact, and **no probe was +1**: 38 of 42 cases were −1 and
4 were 0.

**1-D** (all 26 cases in 29 s; each solve < 0.05 s; exact: DD 9.869604, DN 2.467401):

| BC | form | gen | d | m | λ̂ (cert / re-cert) | bracket | stock |
|---|---|---|---|---|---|---|---|
| DD | L2 | none | 1 / 2 / 3 | 16 / 39 / 72 | uncertain | lo 1.5e−4 / 0.096 / 1.079 | 1.0e−4 / 0.091 / 1.051 |
| DD | L2 | product | 1 | 24 | 1.0000000 (+1/+1) | [0.999991, 1.000047) | 0.99997 |
| DD | L2 | product | 2 | 51 | 9.8692933 (+1/+1) | lo 9.869251 | 9.86899 (stock one-SDP 9.869292) |
| DD | L2 | product | 3 | 88 | 9.8696044 (+1/+1) | lo 9.869553 | 9.86960 |
| DD | L2 | faces | 1 / 2 / 3 | 20 / 45 / 80 | 1.0000000 / 9.8694053 / 9.8696044 (+1/+1 each) | [0.999991, 1.000066) at d = 1 | – |
| DD | energy | product | 1 / 2 | 29 / 58 | 9.7013012 / 9.8696044 (+1/+1) | [9.701248, 9.701850) at d = 1 | 9.70065 / 9.86960 |
| DD | energy | faces | 1 | 24 | **9.7165200** (+1/+1) | [9.716459, 9.717287) | – |
| DD | energy | faces | 2 | 51 | 9.8696044 (+1/+1) | | – |
| DN | L2 | none | 1 / 2 / 3 | 16 / 39 / 72 | uncertain | lo 0 / 0.096 / 0.925 | 1.8e−4 / 0.109 / 1.122 |
| DN | L2 | product | 1 / 2 / 3 | 24 / 51 / 88 | 1.0000000 / 2.4674011 / 2.4674011 (+1/+1) | | 0.99997 / 2.46740 / 2.46740 |
| DN | L2 | faces | 1 / 2 / 3 | 20 / 45 / 80 | 1.0000000 / 2.4674011 / 2.4674011 (+1/+1) | | – |
| DN | energy | product | 1 / 2 | 29 / 58 | 2.4671737 / 2.4674011 (+1/+1) | [2.467162, 2.467219) at d = 1 | 2.46709 / 2.46740 |
| DN | energy | faces | 1 / 2 | 24 / 51 | **2.4673926** / 2.4674011 (+1/+1) | | – |

**2-D, DD × DN, d = 1** (5 min; 0.2–0.7 s per solve; exact: full 12.337006, dir 1 9.869604,
dir 2 2.467401):

| form | gen | sel | m, Gram | λ̂ (cert / re-cert) | bracket | stock |
|---|---|---|---|---|---|---|
| L2 | none | full / 1 / 2 | 512, 1×64 | uncertain (≤ 0.0022) | [0, 0.032) / [0, 0.0072) / [0, 0.0096) | 0.00165 / 1.0e−4 / 3.3e−4 |
| L2 | product | full / 1 / 2 | 1152, 2×64 | uncertain (≤ 0.0036) | [0, 0.090) / [0, 0.039) / [0, 0.018) | 0.00241 / 0.00132 / 0.00048 |
| L2 | faces | full | 768, 5×64 | **2.0000000** (+1/+1) | [1.999944, 2.006721) | 1.99994 |
| L2 | faces | 1 | 768, 5×64 | **1.0000000** (+1/+1) | [0.999972, 1.001779) | 0.99997 |
| L2 | faces | 2 | 768, 5×64 | **1.0000000** (+1/+1) | [0.999972, 1.001629) | 0.99997 |
| energy | faces | full | 920, 5×96 | 12.0412727, not certified (0; +1 at 12.0292314; certified value: lo 12.035914) | [12.035914, 12.156398) | 12.03906 |
| energy | faces | 1 | 844, 5×80 | **9.7165200** (+1/+1) | [9.715480, 9.758855) | 9.71626 |
| energy | faces | 2 | 844, 5×80 | **2.4673926** (+1/+1) | lo 2.467388 | 2.46725 |
| energy | product | full | 1376, 2×96 | uncertain | [0.337, 1.349) | 2.96% (0.365) |

**2-D, d = 2, faces** (6 min; m 3978, 5×324; 22–26 iterations, 35–42 s per solve):

| sel | λ̂ (cert / re-cert) | stock |
|---|---|---|
| full | 12.3368080 (+1 / +1) | 12.33508 |
| 1 | 9.8694055 (+1 / 0) | 9.86930 |
| 2 | 2.4674011 (+1 / +1) | 2.46732 |

Readings (MEASURED unless marked):

- **Shapes.** Every shape equals the stock shape and the reference pricing formulas (c =
  (d+1)(5d+3): none m = 2^(N−1)c^N, faces m = 2^(N−1)[c^N + 2N(d+1)c^(N−1)], Gram
  (2(d+1)²)^N). That includes every energy shape (1-D 24; 2-D 920 and 844).
- **Values.** In the cases whose λ̂ is itself certified (+1: every face-generator row
  except 2-D energy full, and the 1-D product rows), the container reproduces stock's number,
  where it differs it is higher (tight tolerances, objective mode), and it is never above
  exact. The one-SDP value 9.8692933 matches stock's 9.869292. The other rows are not
  covered by that statement:
  - 2-D energy faces full: λ̂ 12.0412727 is not certified; the certified values 12.035914
    (bisection lo) and 12.0292314 (re-certification) are below stock's 12.03906;
  - the no-generator rows and the 2-D product rows are budget-limited brackets, and some
    lie below stock's certified values: 1-D DN L2 none d = 1/2/3 lo 0 / 0.096 / 0.925
    against 1.8e−4 / 0.109 / 1.122; 2-D L2 none full [0, 0.032) against 0.00165; 2-D L2
    product full [0, 0.090) against 0.00241; 2-D energy product full [0.337, 1.349)
    against 0.365.
- **The d = 1 face ceiling.** L2 with faces is 1.0000000 per direction, for DD and DN, in
  1-D and 2-D.
- **Tensor sum rule.** The 2-D full value is 2.0000000 = dir 1 + dir 2. At d = 2 the full
  value 12.3368080 equals 9.8694055 + 2.4674011 to 1.1e−7 relative.
- **B1 = B2.** The 2-D per-direction values equal the 1-D values at the same basis (energy
  DD 9.7165200 in both, to 7e−9 relative).
- **No generators.** Without generators, or with the all-direction product in 2-D, the reach
  is about 0 (triage row C). Those objective optima end uncertain and the brackets are
  budget-limited: soft values, as the reference found.

### 11.4 Measured: 3-D (N = 3, DD × DN × DN, d = 1)

**I0** passes for DD×DN×DN and DN×ND×DD [0 1;1 2;0 3]: 32 checks each, 5.2 s.

**Smoke test** (`'3d0'`): L2 full, no generator, λ = 0, all rows, MOSEK defaults.

- Shape m 16384, 1×512, nnz 2.02M; the build took 5 s.
- Verdict **+1**: rel_b 2.08e−7, psd_relmin +4.5e−11.
- MOSEK 122.2 s (the reference measured 88.6 s), 10 iterations, 24 threads, peak 8.2 GB.

**Pricing** (MEASURED with 2- or 3-iteration capped solves, all rows, tight tolerances):

| family | m, Gram | nnz | build | MOSEK s/iteration | peak | full solve (inferred, 20–27 it) |
|---|---|---|---|---|---|---|
| L2, all 6 faces (full, or per direction: same m and nnz) | 28672, 7×512 | 19.0M | 56 s | ≈ 110 | 13.9 GB | **28–46 min: skipped** |
| energy dir 1, all 6 faces | 31488, 7×640 | 23.4M | 93 s | ≈ 140 | 19.7 GB | **> 45 min: skipped** |
| L2 dir i, faces of direction i only | 20480, 3×512 | 7.7M | 20 s | 25–27 (10 s setup) | 9.7 GB | 10–12 min |
| energy dir i, faces of direction i only | 22528, 3×640 | 9.5M | 32 s | 37–40 (14 s setup) | – | 13–18 min |

The per-direction tests therefore use **the faces of the tested direction only**. The
prediction is unchanged, for two reasons (DERIVED, reference study):

- The B2 certificate, the 1-D certificate ⊗ T_î², needs only g_i.
- B1 holds for any generator set.

That cone is also a **subcone** of the all-faces cone. So the sum of the three certified
per-direction points is a certificate for the full all-faces L2 problem at Σλ_i:
Σ_i Q_i(λ_i) = Q_full(Σλ_i), G0 blocks add, and each face block comes from its own
direction. `'3dsum'` checks this sum on the all-faces SDP with one mat-vec and the
`heatNd_solve` gate, with no solve.

**Lower-side tests, plus one upper side (L2 dir 1 at 1.10)**, of the 1-D-predicted values
(`'3d'`, one λ per run, tight tolerances, all rows). Every 1.01 upper side ended MOSEK
UNKNOWN, which carries no information, so those upper sides are untested. Predictions from
the 1-D objective at the same basis (seconds): L2 DD 1.0000000, L2 DN 1.0000000, energy DD
9.7165200, energy DN 2.4673926. The last row is the final review's own solve
(`…/scratchpad/heatNd/review_poincare/log_3d_e999.txt`), not a `'3d'` run.

| test | λ | verdict | rel_b / psd_relmin | MOSEK | iterations | peak |
|---|---|---|---|---|---|---|
| L2 dir 1 (DD) | 0.99 | **+1** | 1.3e−7 / +2.4e−12 | 607 s | 20 | 9.7 GB |
| L2 dir 1 | 1.01 | 0 (UNKNOWN; rel_b 1.6e−2) | – | 497 s | 17 | 9.7 GB |
| L2 dir 1, MOSEK defaults | 1.01 | 0 (UNKNOWN; the same iterate, rel_b 1.6e−2) | – | 468 s | 17 | 9.7 GB |
| L2 dir 1 | 1.10 | **−1** (verified Farkas ray, relative violation 0) | – | 494 s | 18 | 9.6 GB |
| L2 dir 2 (DN) | 0.99 | **+1** | 1.2e−8 / +3.0e−14 | 628 s | 20 | 9.7 GB |
| L2 dir 3 (DN) | 0.99 | **+1** | 1.4e−8 / +6.1e−14 | 635 s | 20 | 9.7 GB |
| **L2 full, all faces, sum of the three** | **2.97** | **+1** (no solve) | 8.2e−8 / +5.2e−14 | – | – | 4.5 GB |
| energy dir 1 (DD) | 9.6193548 (0.99 ×) | **+1** | 1.7e−8 / +1.5e−13 | 855 s | 19 | 11.2 GB |
| energy dir 1 | 9.8136852 (1.01 ×) | 0 (UNKNOWN; rel_b 1.4) | – | 784 s | 19 | 11.2 GB |
| energy dir 2 (DN) | 2.4427187 (0.99 ×) | 0: MOSEK stopped at the 1080-s cap (UNKNOWN); its iterate has rel_b 5.4e−9, psd_relmin +8.8e−15 | – | 1086 s | 25 | 11.3 GB |
| energy dir 2, MOSEK defaults | 2.4427187 (0.99 ×) | **+1** | 1.5e−7 / +1.1e−12 | 571 s | 13 | 11.2 GB |
| energy dir 1 (DD), MOSEK defaults, final review | 9.7068035 (0.999 ×) | **+1**; semantically checked: Gram PSD, G(d) vs analytic Q 2.0e−5, Ritz min +0.16796 against +0.16829 analytic | 2.04e−7 / +5.5e−12 | 480.9 s | 13 | 11.9 GB |

All solves used 24 MOSEK threads. No other MATLAB process was running at the start of any
3-D run (the runner checks; not monitored throughout).

**Readings.**

- MEASURED: I0, the smoke test, and certified lower sides for every per-direction L2 test
  at 0.99 × prediction. Also certified: the full L2 at 2.97 = 0.99 × 3.0000 in the all-faces
  cone, which is the heat benchmark's generator set.
- The upper side of L2 dir 1 at 1.01 gave no certificate (MOSEK UNKNOWN, with tight and
  with default tolerances) rather than a certified −1, so at 1.01 it is untested. At 1.10
  the verdict is a certified −1, so the 3-D DD L2 per-direction constant is certified in
  [0.99, 1.10) around the predicted 1.0000: resolution 0.99 below, 1.10 above.
- Energy dir 1, the DD direction that sets the heat gap, is certified at 0.99 × e_DD
  (this run) and at 0.999 × e_DD = 9.7068035 (final review). The upper side at 1.01 × e_DD
  is untested (UNKNOWN). B1 with the 1-D certified −1 at 9.7172870 gives the 3-D value
  < 9.71729 (DERIVED, final review: the face weights of the other directions fold into
  the PSD weight). So **3-D e_DD ∈ [9.7068, 9.7173]**, consistent with the 1-D e_DD
  9.7165200 (B1 = B2) at 0.1% resolution.
- Energy dir 2 at 0.99 × e_DN is certified (+1) with MOSEK defaults, in 13 iterations. The
  tight solve stopped at the time cap with an iterate that already met the residual and PSD
  gates; the rule does not coerce that into +1. The upper side for DN (1.01 × e_DN = 2.4921)
  lies above μ_DN = 2.4674, so it is a soundness test, not a sharpness test, and was not run.
- NOT MEASURED, priced over 20 min:
  - the full L2 upper side (3.03);
  - the per-direction tests with all six faces;
  - the full energy form.
- Also not run: L2 dir 2 and 3 at 1.01, and energy dir 3. Dir 3 is the image of dir 2
  under s2 ↔ s3, both DN on [0,1], with an identical basis, so its SDP is a permutation of
  dir 2's (DERIVED, not measured).

### 11.5 Decision table (reference study) and how to use it

Run the steps in order. Each is cheaper than the next, and a failure stops the triage.

1. **I0** (`heatNd_poincare_ops(pie,true)`, seconds at any N). A failure implicates T/A.
2. **References** at the SAME per-direction basis and generators as the LPI under test:
   R_i is the 1-D (or lower-dimensional) per-direction λ̂. The tensor bounds are:
   - B1: λ̂_i(N-D) ≤ R_i;
   - B2: λ̂_i(N-D) ≥ R_i and λ̂_full ≥ Σ_i R_i, when T_j is in the direction-j basis and
     each generator depends on one direction only (faces do; the all-direction product
     does not).

   "Reaches" means within about 2e−4 relative, the gate resolution.
3. **Triage rows**, checked first:
   - **A.** λ̂ > exact (μ_i or λ1), or κ̂ > λ1, beyond 2e−4: UNSOUND. Suspect T/A (a BC
     error changes the spectrum) or a verdict gate that is too loose. The 1.001 × exact
     probe must never be +1.
   - **B.** −1 `deficient` (0 = b ≠ 0 on a row with no variable) at every λ > 0: the basis
     misses a target monomial (tensor d = 0 misses the DD monomial sθ; test (7)). Degree
     or basis, not numerics.
   - **C.** λ̂ ≈ 0 without generators or with the all-direction product, while λ = 0
     certifies: EXPECTED, since there is no interval multiplier. The lever is face
     generators, not a bug.
4. **Decision table.** Assumes I0 passes and generators are present. H = the heat LPI
   reaches its reference; P = the full Poincaré test reaches its reference; Pi = every
   per-direction test reaches its reference.

| row | H | P | Pi | reading |
|---|---|---|---|---|
| 1 | yes | yes | yes | consistent |
| 2 | any | yes | some no | violates B1 (full ≤ Pi + Σ_{j≠i} μ_j < λ1): the per-direction harness is wrong (A_i built differently from A0, or another basis or generator set), or the full test is unsound (row A) |
| 3 | yes | no | yes | the full value is below Σ R_i: the full-target assembly (sum or ctranspose of the containers), or generators missing in some direction; T/A are fine |
| 4 | yes | no | no | the heat P search compensates with P ≠ T. This is normal at low degree or without generators. The Poincaré test is uninformative at that basis: raise d, add faces, or use the energy form |
| 5 | no | yes | yes | if T is in the heat P family, heat R represents (1−ε²)T*T, and heat Q contains the Poincaré basis and generators, then P = T is feasible and heat must certify ≥ P̂ − r. A failure is then a heat-assembly bug ((15a) P*T = T*P, 'symmetric' on a non-self-adjoint expression, the sign or scale of kP*T, ε) or numerics. If the premise fails: the P basis or degree |
| 6 | no | yes | some no | contradicts B1, as row 2 |
| 7 | no | no | yes | the full A0 assembly or cross-direction terms (recheck I0 on A) |
| 8 | no | no | no | (a) I0 fails: T/A; (b) deficient: degree; (c) no generators: Psatz; (d) all present at d ≥ 2: numerics (MOSEK status, rel_b, psd, iterations, raggedness). If the 1-D/2-D reference at the same basis also fails, it is basis or degree, not dimension |

Raggedness means a failure below a certification, or a MOSEK UNKNOWN band just below the
threshold. It is boundary numerics: report the bracket, and do not attribute it to structure.

### 11.6 Applied to the heat benchmark's 3-D result

The heat bench 3-D result (d = 0) is κ̂ ∈ [14.60, 14.70) (§6.3), that is ≥ 0.986 κ\*, with
gap k\* − k̂ ∈ (0.104, 0.204]. Step by step:

1. **I0 passes at N = 3** (MEASURED). T and A, which the heat LPI uses unchanged, are
   correct. Rows A (the T/A half), 7 and 8(a) are excluded.
2. **P and Pi against their references in the L2 form**, at the heat's generator set. In
   1-D and 2-D they reach them at the gate resolution (~2e−4, §11.5): the ceiling 1.0000
   per direction and the 2-D sum rule 2.0000. The 3-D evidence is coarser:
   - per direction, certified ≥ 0.99 × 1.0000; the upper side is untested (1.01: MOSEK
     UNKNOWN), except dir 1, certified < 1.10 (MEASURED, faces of direction i);
   - full, certified ≥ 2.97 in the all-faces cone (MEASURED, summed certificate); upper
     side not run.

   At that resolution nothing contradicts B1 or B2 in the 3-D container machinery both
   LPIs share: `poscopvar` and its face terms (then `heatNd_posw`; bit-identical SDPs, §9), `lpi_eq_cdopvar 'symmetric'`, the 3-D
   cell indexing, and the `heatNd_sdp`/`heatNd_solve` route. Rows 2, 3 and 6 are not
   indicated, reading "reaches" in 3-D as "certified ≥ 0.99", not as ~2e−4.
3. **H.** The heat bench never reaches κ\*, at any N. Its gap is 0.1531 in 1-D and
   0.1530–0.1532 in 2-D (§6). Read against κ\* the 3-D result is "H no, P yes, Pi yes"
   (P, Pi at 0.99 resolution), i.e. row 5.
   - Row 5's premise fails for bench: T ∉ Z_0 (and −A0 ∉ Z_0 for N ≥ 2).
   - So row 5 points at the heat P basis/degree, not at an assembly bug.
   - Its implication is satisfied anyway: heat κ̂ ≥ 14.60 is far above the L2 form's
     P̂ = 3.0 at this basis.
   - The L2 form at d = 1 is weak (the 1.0000 face ceiling). It checks structure, not
     reach.
4. **Energy form**, the informative control:
   - The per-direction energy constants at the heat's generator set and the tensor d = 1
     basis are e_DD = 9.7165200 and e_DN = 2.4673926 (MEASURED, 1-D; 2-D equal to 7e−9).
   - The heat bench κ̂ matches **Σ_i e_i** in every N measured:
     - 1-D: 9.7165200 lies in the heat bracket [9.716459, 9.716534);
     - 2-D: 9.7165200 + 2.4673926 = 12.1839126 lies in [12.183790, 12.183978), and in
       the union bracket [12.183809, 12.183978) of §6.2;
     - 3-D: Σ = 14.6513052 lies in [14.60, 14.70). This is the value extrapolated from
       the lower dimensions, the fix report's gap prediction 14.651 reached by another
       route: heat = Σ e_i was measured in 1-D and 2-D, so Σ of the 1-D e_i equals λ1
       minus the 1-D/2-D heat gap by construction. It adds no independent evidence about
       3-D, whichever LPI supplies it.
   - So the heat bench gap, 0.1531, numerically equals the DD direction's energy-form loss
     at this basis, μ_DD − e_DD = 0.1531. The DN directions lose ≤ 1e−5 each.
   - In 3-D the DD energy constant is certified ≥ 0.99 e_DD = 9.6193548 (MEASURED, this
     study) and ≥ 0.999 e_DD = 9.7068035 (MEASURED, final review, §11.4). Its upper side
     is untested (1.01: UNKNOWN); B1 bounds it by 9.71729 (DERIVED). So **3-D e_DD ∈
     [9.7068, 9.7173]**, and a 3-D-specific DD energy loss is ≤ 0.0097, 6% of the heat
     gap 0.153. At this study's own 0.99 resolution it could have been up to 0.097 (63% of
     the gap).
   - The bound on a 3-D-specific component of the **heat** shortfall comes from the heat's
     own bracket, not from the Poincaré tests: the 3-D gap (0.104, 0.204] differs from the
     extrapolated 0.1531 by (−0.049, +0.051], i.e. ±0.05, up to 1/3 of the gap.
   - INFERRED: the mechanism is not derived. In 1-D, −A0 = −I is in Z_0, and the bench
     Q basis (`Tspan` dp = 0, which is tensor d = 1 in DD) contains the energy basis, so
     heat ≥ e_DD there is plausible. For N ≥ 2 the additivity is an observation, to ≤ 1e−4.
5. **Conclusion, with its resolution.** T/A are right (I0). The shared container machinery
   shows nothing in 3-D that contradicts the tensor bounds, at 0.99 resolution. The 3-D DD
   energy constant equals the 1-D e_DD to 0.1% (final review), so there is no 3-D-specific
   loss in it larger than 0.0097. The 3-D heat bracket contains the value extrapolated from
   the lower dimensions, and it bounds any 3-D-specific component of the heat shortfall by
   ±0.05 (up to 1/3 of the gap). Nothing tighter is measured.

**What the diagnostic would flag if the heat LPI regressed** (re-run `test_heatNd_poincare`
and `heatNd_poincare_validate('I0')` first, in seconds):

- **The heat κ̂ exceeds λ1**, or a 1.001 × exact Poincaré probe certifies: row A, UNSOUND.
  Look at the verdict gate or at T/A.
- **I0 fails**: T/A (`heatNd_pie`), the BC map or the cell indexing. Stop; neither LPI can
  be trusted.
- **I0 passes, but a Poincaré per-direction value moves off its 1-D reference**, or the
  2-D/3-D sum rule fails: the shared container machinery (rows 2/3/6), i.e. `poscopvar`,
  its face codes, `lpi_eq_cdopvar 'symmetric'`, sorted-variable cells, or `heatNd_sdp`. Fix
  that before touching the heat LPI.
- **I0 and all Poincaré references pass, but heat κ̂ drops below Σ_i e_i** by more than its
  bracket (for example 1-D below 9.71646, 2-D below 12.1838, or 3-D certified −1 below
  14.65): suspect a heat-specific regression (row 5/8(d)). Σ_i e_i is an empirical
  reference, not a bound, so a drop is a flag to investigate, not a proof. Candidates:
  - the P basis (Z_d), the (15a)/(15b) encoding ('split', 'symmetric' pairing), the sign or
    scale of kappa·X2, ε;
  - the `Tspan` R/Q bases, or the faces on R and Q;
  - numerics: rel_b at the gate, UNKNOWN bands, the row set. Retry with `tight`/`rows`.

  If this happens only in 3-D while 1-D/2-D still match, the 3-D Poincaré tests exonerate
  the shared 3-D machinery (to their resolution: 0.99, and 0.999 for DD energy) and point
  at the heat's own 3-D build or MOSEK accuracy (3-D rel_b floor 3.3e−7, §6.3).
- **Heat κ̂ rises above Σ_i e_i**: not a failure, but it breaks the empirical additivity
  above. Re-examine the P search before reading it as progress.

### 11.7 Measured vs inferred, and open items

MEASURED:

- Every I0 error, shape, rank, timing, memory peak, verdict, rel_b and λ̂ above.
- The MOSEK thread count (24 every solve).
- The QR blow-up (> 70 GB).
- The per-iteration costs of the all-faces and restricted 3-D families. This corrects the
  reference's inferred ≈ 54 s/iteration (11–29 min per solve) for 3-D faces L2 to
  ≈ 110 s/iteration measured.
- The smoke test took 122 s here, against the reference's 88.6 s.
- The summed 3-D certificate.
- The equality of container and stock shapes, and the container values against stock.

DERIVED: the exact constants, the identities, and the P = T and P = −A0 links to Cor. 35.
The tensor bounds B1 and B2 come from the reference study. The subcone argument for
`'3dsum'`, the invariance of the per-direction prediction under restricting the faces, and
the dir 2 ↔ dir 3 symmetry are derived here.

INFERRED:

- The full-solve prices from the capped solves (20–27 iterations, as in 1-D/2-D).
- That the 3-D L2 per-direction value is exactly 1.0000. Certified: [0.99, 1.10) for
  dir 1, ≥ 0.99 for dirs 2 and 3; the upper side at 1.01 is untested (UNKNOWN).
- That the 3-D energy constants equal the 1-D ones. Certified: ≥ 0.99 × (dir 2); dir 1
  ≥ 0.999 × (final review) and < 9.71729 by B1 (derived), i.e. [9.7068, 9.7173]; the
  upper side at 1.01 × is untested (UNKNOWN).
- The mechanism behind heat κ̂ = Σ_i e_i. It is a measured coincidence to ≤ 1e−4 in 1-D and
  2-D, and consistent in 3-D to ±0.05 (the heat bracket); it is not derived.

Open:

- The 3-D upper sides at 1.01 are untested (MOSEK UNKNOWN, with tight and with default
  tolerances), not certified −1. L2 dir 1 is certified −1 at 1.10.
- The all-faces per-direction, full-upper and full-energy 3-D solves are over the 20-min
  budget. They need a machine or window where 30–50 min solves are allowed, or a smaller
  (`Tspan`-like) Gram basis.
- `heatNd_lindep` cannot be used on the 3-D faces families (memory).
- Both maintainer questions of §10 are answered. The paper's Table 1 d = 0/1 rows, from
  an earlier PIETOOLS version, are not reproduced here either.

## 12. Degree tailoring of R and Q (MMP, 10/08/2026)

Question: does sizing R and Q by the lift/weight rules of the proof program (sopvar_lift_notes
Sec. 7: lower-kernel support {min(i,j) <= D, max(i,j) <= D + 2w + 1}, multiplier degree <= 2w,
lift D and weight w separate, the product Psatz term one degree below the plain term) reach the
same certified rate as the stock presets with a smaller SDP? `heatNd_tailor(N,dlist,bopts,bc)`
answers it: for every P degree d it bisects `bench` and `heavy`, reads Dmin, Dmax and the
multiplier degree off the targets E1 and X1 + X2, and bisects the grid D = Dmin + {0,1},
w = rule + {0,1}, Psatz none / product at w-1 / product at w / faces at w / faces at w-1, each
R and Q from its own target, no joint cap. It uses two new `heatNd_lpi` options: `RQ 'custom'`
(`Rdeg`, `Qdeg` copquadvar specs used as given) and `psatz_offset` (the Psatz terms at `int`
reduced by it; with `product` and 1 this is the Markov-Lukacs pair S0 + g S1, deg S1 = deg S0 - 2).
Defaults reproduce every earlier program. Measured 1-D, eps = 0.1, certified bisection, the
decisive rows re-run at rtol 1e-6 (brackets end "uncertain" only because MOSEK returns UNKNOWN
above lambda1, as in Sec. 6); nx = SDP variables, m = equality rows.

**DD, kappa* = 9.869604.** `heavy` (int = mult = 2, joint cap 3, two face terms at the same degree)
certifies 9.8695954 at nx 1854 / 1859 / 1866 and m 87 / 91 / 96 for d = 0 / 1 / 2. The same rate,
to the last digit of the bisection, comes from R at (D,w) = (2,1), Q at (2,2) and the product term
at w - 1: nx 820 / 825 / 832, m 78 / 82 / 87, i.e. 2.25 times fewer variables and 10% fewer rows at
every d. Smaller variants fall short: R(2,1) Q(2,1) gives 9.8048 with the product pair and 9.8690
(gap 6e-4, nx 1028) with two faces; R(2,1) Q(1,2) gives 9.8620; R(1,2) Q(2,2) gives 9.8695860
(nx 853). Q is the binding operator and needs weight 2 with lift 2; R needs weight 1. The
necessary conditions place D at Dmin = 1 and w at 1 for d = 1; the certificate needs D = 2 and,
for Q, w = 2. The two-face form at the same degrees costs 1763 variables against 825 for the
product pair. `bench` (9.7165 at d = 0, 1) is not improved at equal rate: the smallest tailored
SDP reaching it is R(2,1) Q(2,1) with the product pair at d = 0 (9.8043, nx 428 against 495, but
m 64 against 51).

**DN, kappa* = 2.467401.** `bench` already certifies 2.467346 at nx 303 / 308 / 315 (its Tspan
cap a + c <= 1 keeps 3 of the 4 monomials), and `heavy` the same rate at 1110 to 1122. The
tailored R = Q = (1,1) with the product pair reaches 2.467195 (gap 2.1e-4) at nx 213, m 44 at
d = 1; nothing smaller than `bench` reaches 2.467346.

Without any Psatz term no configuration certifies a rate (as Sec. 6 for `listing`); a product
term at w = 0 does not exist and a face term at w - 1 = 0 certifies nothing in DD.

**2-D, DD x DN, kappa* = 12.337006, d = 0 (10/08/2026, rtol 1e-4, tmax 150-200 s per bisection).**
`bench` certifies 12.183764 (gap 0.153) at nx 40409, m 1242; `heavy` 12.336731 (gap 2.7e-4) at
nx 565209, m 4282. The 1-D pattern does not carry over as such: with Q at (D,w) = (2,2) per
direction the product pair certifies nothing (R at w = 1) or 1.5 (R at w = 2), as Sec. 3 found
for the product generator; the faces are needed. Uncapped faces at w - 1 with Q at weight 2 reach
12.335223 (gap 1.8e-3) at nx 339209 to 360009 whatever R is (R at (1,1) weight 1 suffices), and
uncapped faces at w with Q at weight 2 give nx 1036809, beyond the budget. What does carry over is
the split between the two operators: R at the `bench` degrees (Tspan dp 0) with Q at the `heavy`
degrees (Tspan dp 1), faces at the same degree, certifies heavy's 12.336731 at nx 392409, m 2884,
i.e. 31% fewer SDP variables and 33% fewer rows than `heavy`; the reverse split gives bench's
12.183 (Q binds, R does not), and the faces on Q only certify nothing (as Sec. 3). So in 2-D the
stock graded basis (the Tspan subset caps) is what keeps the degree-2 face terms affordable, and
the saving comes from giving Q and R different degrees, not from the Markov-Lukacs pair.
`heatNd_tailor` runs the grid per direction (support_degrees, rule per direction); the split
configurations were run directly with `RQ 'custom'` and the Tspan specs of `RQ_deg`.

**The tensor term set against the faces (10/08/2026, later).** `heatNd_lpi` psatz `'tensor'`: the
Markov-Lukacs tensor set, each box quadratic `(th_i - a_i)(b_i - th_i)` at `int - offset` in its own
direction and their product at `int - offset` in every direction (`poscopvar_direct` codes `2N+2+i`
and 1). 2-D DD x DN, d = 0, R (1,1), Q (2,2) uncapped, certified bisection with a 200 to 400 s budget:

| terms | kappa certified | gap to kappa* | nx | m | solves |
|---|---|---|---|---|---|
| faces at w - 1 | 12.336731 | 2.7e-4 | 339209 | 4002 | 16, 192 s |
| tensor at w - 1 | 12.326178 | 1.1e-2 | 395785 | 4002 | 15, 390 s |
| tensor at w | 6.17 after 2 solves | | 762889 | 7126 | budget |

The faces at w - 1 reach heavy's 12.336731 (the earlier record of 12.335223 for this configuration
was the smaller solve budget) at 339209 variables, below the split configuration (392409) and
`heavy` (565209); the tensor set is larger and certifies less. The same holds on the 2-D coercive
H-infinity slack (`lpi_programming_sopvar/README.md`, `lpi_ineq_sop`). The N-D default of
`get_lift_degs` is now faces at w - 1 whenever every weight is at least 2.
