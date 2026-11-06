# Class-operation baseline: opvar vs sopvar, dopvar vs sdopvar

Measured 11/05/2026 on `ndopvar` at `edbd7341`, with no local changes under `opvar/`, `sopvar/`,
`SOSTOOLS400/dpvar/` or `converters/operator_converters/`. MATLAB R2025b, Windows 11, Intel
i9-14900KF, 63.7 GB. The sweeps ran in about 10 min.

Files: `private/bench_class_ops.m` (benchmark), `private/analyze_class_ops.m` (summary),
`results.csv` (every timed row), `tables.md` (all summary tables). The two functions sit in
`private/` to stay off the MATLAB path. The benchmark needs `opvar/`, so `results.csv` is the
record once the opvar classes are retired. `scaling/` holds scaling-plot data and figures (operation
time against q, stability-LPI parsing time against the Lyapunov degree), with its own README.

## Scope

**Operators.** Single-block operators on L2^n[0,1],

    (P x)(s) = R0(s) x(s) + int_0^s R1(s,t) x(t) dt + int_s^1 R2(s,t) x(t) dt,

with R0, R1, R2 polynomial of degree d in s and in t. `opvar` and `dopvar` hold P as a 3-PI
operator (dim [0 0; n n]); `sopvar` and `sdopvar` hold it as one block on L2[s]. In the decision
classes every coefficient of R0, R1, R2 is a constant plus a linear combination of decision
variables; q is the number of decision variables and nnz the number of nonzero decision
coefficients.

**Operations.**

| name | expression | classes |
|---|---|---|
| plus, minus | A + B, A − B | all four |
| uminus, scalar | −A, 2.5·A | all four |
| mtimes | A·B, composition | opvar, sopvar |
| mtimes_DA, mtimes_AD | D·A, A·D, decision operator D with constant operator A | dopvar, sdopvar |
| plus_const | D + A | dopvar, sdopvar |
| ctranspose | A', the L2 adjoint | all four |
| horzcat, vertcat, blkdiag | [A B], [A; B], blkdiag(A,B) | all four |
| subsref | A(1:⌊n/2⌋,:), A(1,:) at n = 1 | all four |
| inv | inverse of I + 0.1 K | opvar, sopvar |

Not compared: `diff` and `subs`, whose definitions differ between the classes (the opvar versions
return operators on [x; ∂_s x] and [x; x(s0)], the sopvar versions operators on x), and the
composition of two decision operators, which both classes refuse.

**Equivalence.** Each operator is built once in the legacy class – by `rand_opvar`, or as a dopvar
assembled in the dpvar coefficient layout – and converted once by `opvar2sopvar` or
`dopvar2sdopvar`. `bench_class_ops('check')` applies every operation in both classes and applies
both results to 3 random quadratic inputs at 60 random points, each class from its own
definition: multipoly integrals for opvar, `apply_sopvar` for sopvar. A decision result is first
evaluated at random decision values, by `pi_dpvar_eval` for dopvar and by C = A + Bᵀd for sdopvar.
Over the 22 operations the largest relative difference is 3.7e-16.

## Method

- **Time:** one warm-up call, then the median of 3–25 calls filling at least 0.3 s.
- **Storage:** `whos` of the result. The `uminus` result has the size of its operand. `whos`
  counts an array shared copy-on-write once per holder.
- **Exponents:** over (n, d), least squares log t = c + α log n + β log d. Over q, the slope of
  log t between the two largest q.

| part | classes | grid | decision structure |
|---|---|---|---|
| const | opvar, sopvar | n ∈ {1,2,4,8} × d ∈ {1,2,4,8} | none |
| dvarB | dopvar, sdopvar | n ∈ {1,2,4,8} × d ∈ {1,2,4} | each coefficient its own variable, q = n²((d+1) + 2(d+1)²) ≤ 3520 |
| dvarQ | dopvar, sdopvar | n = d = 2, q ∈ {10², …, 10⁶} | each coefficient one variable drawn from q; nnz = 168 for every q |
| dvarS | dopvar, sdopvar | n = d = 2, q ∈ {10³, …, 10⁶} | as dvarQ, both operands on one list |
| dvarN | dopvar, sdopvar | n = d = 2, q ∈ {10², …, 10⁵} | each coefficient ⌈4q/84⌉ variables drawn from q; nnz ≈ 3.9q |
| inv | opvar, sopvar | n ∈ {1,2,4} × d ∈ {1,2} | none |

In dvarB, dvarQ and dvarN the second operand of a binary operation has its own decision variables;
the two lists are disjoint.

## Results

### Constant operators, n = d = 8 (const)

| operation | opvar (s) | sopvar (s) | sopvar / opvar | ratio over the grid [min, max] | opvar α, β | sopvar α, β |
|---|---|---|---|---|---|---|
| plus | 1.61e-3 | 5.05e-4 | 0.31 | [0.24, 0.58] | 0.07, 0.09 | −0.02, 0.05 |
| minus | 1.93e-3 | 4.69e-4 | 0.24 | [0.24, 0.45] | 0.05, 0.09 | −0.03, 0.01 |
| uminus | 1.45e-4 | 1.26e-5 | 0.09 | [0.04, 0.09] | −0.02, 0.00 | 0.08, 0.02 |
| scalar | 2.59e-4 | 2.48e-5 | 0.10 | [0.07, 0.10] | −0.02, −0.01 | 0.04, 0.04 |
| mtimes | 0.252 | 0.0266 | 0.11 | [0.03, 0.11] | 0.30, 0.73 | 0.53, 0.52 |
| ctranspose | 4.99e-4 | 1.08e-3 | 2.16 | [0.50, 2.16] | −0.02, −0.07 | 0.17, 0.15 |
| horzcat | 1.72e-3 | 4.80e-4 | 0.28 | [0.28, 0.42] | 0.10, 0.11 | 0.01, 0.03 |
| vertcat | 1.66e-3 | 4.64e-4 | 0.28 | [0.27, 0.42] | 0.11, 0.10 | 0.00, 0.02 |
| blkdiag | 5.03e-3 | 5.27e-4 | 0.10 | [0.10, 0.12] | 0.06, 0.09 | 0.05, 0.06 |
| subsref | 9.79e-4 | 2.40e-4 | 0.25 | [0.25, 0.39] | 0.07, 0.06 | −0.02, 0.00 |

### Decision operators: one decision operand

q = 10⁶ with nnz = 168 (dvarQ), and q = 10⁵ with nnz = 3.9·10⁵ (dvarN). Times in seconds.

| operation | dvarQ dopvar | dvarQ sdopvar | q-exponent dopvar, sdopvar | dvarN dopvar | dvarN sdopvar | q-exponent dopvar, sdopvar |
|---|---|---|---|---|---|---|
| uminus | 2.93e-4 | 1.40e-5 | −0.02, 0.02 | 1.53e-3 | 1.07e-3 | 0.66, 1.58 |
| scalar | 1.97e-3 | 1.52e-5 | 0.00, 0.03 | 2.69e-2 | 9.99e-4 | 0.82, 1.36 |
| plus_const | 2.63 | 4.23e-4 | 1.19, 0.02 | 0.182 | 4.06e-4 | 1.01, 0.01 |
| mtimes_DA | 3.02 | 3.05e-3 | 1.04, 0.00 | 2.90 | 0.532 | 0.92, 0.93 |
| mtimes_AD | 2.94 | 4.49e-3 | 1.02, 0.02 | 3.40 | 0.589 | 0.93, 0.91 |
| ctranspose | 7.65e-4 | 3.65e-4 | 0.00, −0.01 | 1.07e-2 | 3.35e-3 | 0.80, 0.84 |
| subsref | 2.85e-2 | 1.65e-4 | 0.84, −0.01 | 2.75e-2 | 3.35e-4 | 0.90, 0.24 |

### Decision operators: two decision operands

q = 10⁶ with nnz = 168: disjoint lists (dvarQ) and one shared list (dvarS). Times in seconds.

| operation | disjoint: dopvar | disjoint: sdopvar | sdopvar / dopvar | q-exponent dopvar, sdopvar | shared: dopvar | shared: sdopvar | sdopvar / dopvar | q-exponent dopvar, sdopvar |
|---|---|---|---|---|---|---|---|---|
| plus | 1.58 | 3.28 | 2.07 | 1.06, 1.20 | 1.76 | 0.0536 | 0.03 | 1.10, 0.99 |
| minus | 1.62 | 3.19 | 1.97 | 1.06, 1.23 | 1.63 | 0.0489 | 0.03 | 1.06, 0.90 |
| horzcat | 3.21 | 3.28 | 1.02 | 1.13, 1.22 | 3.31 | 0.0488 | 0.01 | 1.15, 0.96 |
| vertcat | 3.21 | 3.44 | 1.07 | 1.18, 1.24 | 3.37 | 0.0482 | 0.01 | 1.18, 0.96 |
| blkdiag | 3.48 | 3.22 | 0.92 | 1.15, 1.18 | 3.52 | 0.0776 | 0.02 | 1.18, 0.90 |

With nnz ≈ 3.9q (dvarN, disjoint lists, q = 10⁵) the sdopvar/dopvar ratios are 0.57 (plus), 0.69
(minus), 0.74 (horzcat), 0.82 (vertcat) and 0.70 (blkdiag), with q-exponents 1.14–1.42 in both
classes. For the free-variable structure (dvarB) at n = 8, d = 4 (q = 3520) the ratios are
0.39–1.08 for the two-operand operations and 0.005–0.62 for the others. The sdopvar exponents α of
the two-operand operations (0.49–0.73) exceed the dopvar ones (0.14–0.29), consistent with
finding 5: here q grows as n² and the lists are disjoint.

### Inverse (inv)

| n | d | opvar (s) | sopvar (s) | sopvar / opvar | residual opvar | residual sopvar |
|---|---|---|---|---|---|---|
| 1 | 1 | 0.0992 | 0.0166 | 0.17 | 4.0e-10 | 2.7e-11 |
| 1 | 2 | 0.109 | 0.0255 | 0.23 | 6.0e-10 | 5.6e-11 |
| 2 | 1 | 0.113 | 0.0174 | 0.15 | 6.9e-10 | 7.1e-11 |
| 2 | 2 | 0.153 | 0.0267 | 0.17 | 7.7e-10 | 1.0e-10 |
| 4 | 1 | 0.204 | 0.0238 | 0.12 | 5.5e-10 | 1.7e-10 |
| 4 | 2 | 0.320 | 0.0584 | 0.18 | 1.5e-9 | 3.2e-12 |

Residual: sup over 3 random quadratic inputs f and 60 points of |P (P⁻¹ f) − f| / sup |f|.

### Storage (bytes of the operator)

| part | range | new / legacy |
|---|---|---|
| const | n ≤ 8, d ≤ 8 | 0.45–0.93 |
| dvarB | n ≤ 8, d ≤ 4 | 0.41–0.64 |
| dvarQ | q = 10²–10⁶, nnz = 168 | 0.33–0.43 |
| dvarN | q = 10²–10⁵ | 0.42–0.49 |

At q = 10⁶ with nnz = 168 the sdopvar takes 133.8 MB and the dopvar 401.3 MB. The coefficients
occupy a few kB of this; the rest is the decision variable name list, 134 bytes per name for the
names 'a1' … 'a1000000'. The
sdopvar holds one list; each of the three kernel dpvars of the dopvar holds one, and `whos`
counts a list shared among them once per kernel, so 401.3 MB bounds the dopvar from above.

## Findings

1. **Constant operators.** Over the whole grid the sopvar time is 0.03–0.58 of the opvar time for
   every operation except the adjoint. Composition shows the largest absolute gap: at n = d = 8,
   0.0266 s against 0.252 s. The adjoint ratio runs from 0.50 (n = 2, d = 1) to 2.16
   (n = d = 8): the sopvar adjoint grows with size (α = 0.17, β = 0.15), the opvar adjoint does
   not (α = −0.02, β = −0.07).
2. **Overhead regime.** For n ≤ 8 and d ≤ 8 every constant-operator operation except composition
   has |α|, |β| ≤ 0.17 in both classes. Those times measure the fixed cost per call, not the
   arithmetic, and their exponents do not extrapolate to larger n or d. Composition grows with
   α = 0.30, β = 0.73 (opvar) and α = 0.53, β = 0.52 (sopvar).
3. **One decision operand, nnz fixed.** The sdopvar times of uminus, scalar, plus_const, both
   compositions with a constant operator, the adjoint and indexing do not grow with q (exponents
   −0.02 to 0.03 up to q = 10⁶). The dopvar times of plus_const (exponent 1.19), the compositions
   (1.02–1.04) and indexing (0.84) grow with q. At q = 10⁶ the compositions take 2.9–3.0 s in
   dopvar and 3.1–4.5 ms in sdopvar.
4. **One decision operand, nnz ∝ q.** Composition is linear in nnz in both classes (exponents
   0.91–0.93); at nnz = 3.9·10⁵ the sdopvar takes 0.17–0.18 of the dopvar time (0.53 s against
   2.90 s for D·A, 0.59 s against 3.40 s for A·D).
5. **Two decision operands.** With disjoint lists both classes merge the two decision variable
   lists, with q-exponents 1.06–1.24. At q = 10⁶, sdopvar plus and minus take 2.0–2.1 times the
   dopvar time and the concatenations 0.92–1.07 times. With one shared list, sdopvar
   takes 0.01–0.03 of the dopvar time (49–78 ms against 1.6–3.5 s), with q-exponent 0.90–0.99.
   The code accounts for the difference as follows (read from the source, not profiled):
   `@sdopvar/private/sync_basis` compares the lists with `isequal` (line 134) and stops there
   when they are equal; otherwise it forms `unique` of their concatenation (line 140) and
   re-indexes each operand by `ChangeDecVar`, which calls `ismember`. `@dpvar/plus` concatenates
   the two lists (line 184) and merges them with `combine` (line 251), equal or not.
6. **Inverse.** For n ≤ 4, d ≤ 2 the sopvar inverse takes 0.12–0.23 of the opvar time, with
   residual at most 1.7e-10 against at most 1.5e-9.
7. **Storage.** sopvar takes 0.45–0.93 of the opvar storage and sdopvar 0.33–0.64 of the dopvar
   storage. At large q the decision variable names dominate both: 134 bytes per variable per list.

## Limits

- **One spatial dimension.** opvar and dopvar are 1-D classes, so the spatial-dimension axis of
  the sopvar classes has no counterpart here; opvar2d is a separate class.
- **Single blocks.** Operators with R^m components are containers (`copvar`, `cdopvar`);
  `sopvar/Testfolder/Bench_copvar_vs_legacy.m` (09/18/2026) compares those with opvar for plus,
  adjoint and composition only.
- **Synthetic structure.** Kernels are dense random polynomials. Decision structures are synthetic
  (one, or ⌈4q/84⌉, random variables per coefficient); `lpivar` and `poslpivar` produce other
  sparsity patterns.
- **What the numbers cover.** Times are warm medians on one machine; first-call JIT costs are
  excluded. Storage is `whos` of the result, not peak memory.
- **Range.** n ≤ 8, d ≤ 8 for constant operators; q ≤ 10⁶ for decision operators at n = d = 2.

## Reproduce

Run `pietools_path_update` from the PIETOOLS root, then make
`baselines/2026-11-05_class_ops/private` the current folder (the benchmark also uses
`rand_opvar`, `apply_sopvar` and `pi_dpvar_eval` from the path):

    bench_class_ops('check')
    Rc = bench_class_ops('const',[1 2 4 8],[1 2 4 8]);
    Rb = bench_class_ops('dvarB',[1 2 4 8],[1 2 4]);
    Rq = bench_class_ops('dvarQ',[1e2 1e3 1e4 1e5 1e6],2,2);
    Rs = bench_class_ops('dvarS',[1e3 1e4 1e5 1e6],2,2);
    Rn = bench_class_ops('dvarN',[1e2 1e3 1e4 1e5],2,2);
    Ri = bench_class_ops('inv',[1 2 4],[1 2]);
    save('full.mat','Rc','Rb','Rq','Rs','Rn','Ri');  analyze_class_ops('full.mat','.')
