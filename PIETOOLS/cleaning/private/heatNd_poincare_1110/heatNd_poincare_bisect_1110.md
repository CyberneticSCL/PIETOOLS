# test_heatNd_poincare: first bad commit e3a1f667 (bisection, 11/10/2026)

Copied here on 11/19/2026 from the side-session notes of local_485f38f6, with the two scripts it
names (maintainer: "send any unlogged conclusions to repo cleaning").

Side session started by local_485f38f6. The bisection below was report only. On the maintainer's
answer "Yes, do both", the library fix and the test-side count change were then made in the
working tree (section "Fix", next), committed 11/14/2026 as 590876d7 (maintainer: "Yes, go ahead")
and pushed to origin/ndopvar the same day.

## Fix (11/10/2026; commit 590876d7, 11/14/2026)

**Files changed in the working tree** (five code files and one README entry):

- `lpi_programming_sopvar/misc/private/collect_eq_rows.m`. In `drop_multiples`, one line is added
  after `act = (1:nt)';`:

  ```matlab
  act = act(~(cnt==1 & fr==1));
  ```

  - A column whose only compared entry is row 1 – an equation without a decision variable – now
    takes no part in the grouping. Every such equation is written, as before e3a1f667.
  - The key begins with (count, first row), so a group that holds such a column holds only such
    columns. The change therefore removes whole groups; every other group, representative and
    dropped equation stays as it was (from the code).
  - Also changed: the help (description and `drop_multiples`) and a header-log entry.
- `lpi_programming_sopvar/misc/lpi_eq_cdopvar.m` and `lpi_eq_sdopvar.m`: help and log only.
- `lpi_programming_sopvar/README.md`: a log entry with the measurements below. The same file
  holds the parent session's uncommitted hunk (the cx_exec paths under "Tests"), so stage the two
  hunks separately.
- `testing_codes/private/PIETOOLS_demos/sopvar_demos/heatNd_poincare.m`: `build_family` turns the
  row record on (`prog.sopeq = {}`) and returns `R.shape.ndup`, the number of equations left out
  as multiples. Help (OUTPUT `shape`, Cost) and log updated.
- `testing_codes/private/PIETOOLS_demos/sopvar_demos/test_heatNd_poincare.m`: section (3) compares
  `m + ndup` with the pricing formulas and the stock counts. The old lines are kept as (was) lines.

**Measured.** All runs use `git archive` exports of 8c416566: HEAD, and HEAD plus the five code
files. Each run starts on a clean path with a check that the files resolve inside the export, in
one MATLAB process, serially.

| check | HEAD | patched |
|---|---|---|
| `test_heatNd_poincare` | error, affine 6.7e-01 | pass, 44 of 44 checks, 29.3 s |
| (7) family: m, rows without a variable, b third-point error | 6, 1, 6.7e-01 | 10, 5, 5.9e-17 |
| (7) family against the e56c7981 build (At, raw b, three lam) | m differs | bitwise equal |
| (3): m written + ndup = formula or stock count, 2-D | – | 481 + 31 = 512, …, 824 + 20 = 844 |
| stability LPI of E4 (`heatNd_lpi`, 'bench', d = 0), m at N = 1, 2, 3 | 47, 1,181, 29,774 | the same |
| the same LPIs: ndup; rows without a variable | 4, 61, 1,996; 0 | the same |
| the same LPIs: m, nnz and five weighted sums of At and b | – | equal to 17 digits |
| 3-D (q = 2.7M): writer stages eq15a + eq15b, build, peak private | 2.12 + 2.45 s, 31.5 s, 4,221 MB | 1.78 + 2.22 s, 27.0 s, 4,212 MB |
| `bench_lpi_eq_sop`, `lpi_eq_cdopvar` at q ≈ 1e6, N = 1, 2, 3 | 1.070, 1.109, 1.134 s | 1.054, 1.072, 1.109 s |
| the same, peak memory (profiler) and rows | 45.8, 45.9, 47.0 MB | the same |

On the patched export these also pass:

- `test_lpi_eq_sop` (29 checks);
- `test_endtoend_sop` (E1, E2);
- `test_lpi_ineq_sop` (6);
- `test_copquadvar_faces` (97);
- `test_heatNd_lpi` (74);
- `test_executives_sop(false)` (20, including the 2-D stability executive at nx 103,744).

**Cost.** The added line is a mask on two per-equation arrays, O(number of equations).

- In time and peak memory it is within the run-to-run spread at q = 1e4 to 1e6, N = 1 to 3
  (benchmark), and at q = 2.7M in 3-D (measured).
- Rows written change only in programs that have equations without a decision variable. Each such
  equation has b_k ≠ 0 and every coefficient of d at most 5e-13 |b_k|. With exact zeros it reads
  0 = b_k, so the program is infeasible unless b_k is round-off (inferred). There every one of
  them is written, as before 11/01.

**Not run.**

- E4 (15–20 min): its program is unchanged (equal m and fingerprints above), so its result cannot
  change (inferred).
- E3: `e2e_stability_2d` is not in `tests/`.
- No run in the working tree: the parent session was active, with a MATLAB process running. The
  five code files in the working tree are byte-identical to the tested export.

**Open (11/19/2026).**

- `testing_codes/README.md`: the "cause not traced" sentence (c511243a) now names e3a1f667 and
  590876d7 (11/19/2026, in the working tree, uncommitted).
- `test_heatNd_poincare` last ran on 11/14/2026, on an export of 7670fb57 with the fix: 44 of 44
  (message of 590876d7). Later commits changed `poscopvar` and `lpivar_cdopvar` (7ff7e201,
  faca4e51). Their notes say the original call forms, which `heatNd_poincare` uses, are kept, so
  the test is inferred unaffected; not measured.

## Result

- **First bad commit: e3a1f667** (11/04/2026, "LPI equalities: do not write rows that are multiples
  of a written row"). Its only parent, e56c7981, passes all 44 checks. e3a1f667 stops with
  `heatNd_poincare:affine`, third-point error 6.7e-01 (measured).
- **Mechanism.** `drop_multiples` was added by that commit as a local function of
  `lpi_programming_sopvar/misc/private/collect_eq_rows.m`.
  - It groups equations by their whole coefficient column, including the constant term (the
    right-hand side).
  - Of each group it writes only the member with the largest scale.
  - Every equation without a decision variable, 0 = b_k, is a multiple of every other such equation.
    So only the one with the largest |b_k(λ)| is written, and which one that is changes with λ.
- **Where: library code.** The test does not assume a b scaling.
  - Raw b = `D.b*D.bscl` is the concatenated `prog.expr.b`. For the same family at e56c7981 it is
    affine to 5.9e-17.
  - The family construction in `heatNd_poincare` assumes the written rows do not depend on λ.
    e3a1f667 broke that for equations without a decision variable, and the check at line 295
    detects it correctly.
  - For one program the output is still valid: the dropped equations are verified multiples. One
    written 0 = b_i with b_i ≠ 0 keeps the program infeasible, as the five did (inferred; the (7)
    verdict itself was not measured past the guard).
- **Masked failures.** The same commit also fails eight shape checks of section (3): in 2-D, the
  written row counts fall below the pricing formulas. The (7) error hides them. With (7) caught,
  35 of 44 checks pass at e3a1f667 and at 864a00d8 (measured).

## Runs

Each run starts on a clean path:

- `restoredefaultpath`, then the export's `pietools_path_update`, then `testing_codes/private` when
  it exists.
- A guard that the test, the `heatNd_*` files, `lpi_eq_cdopvar`, `poscopvar`, `copquadvar`,
  `lpiprogram` and `sossolve` resolve inside the export, plus `heatNd_path`'s own assertion.
- One MATLAB process per export, run serially.

Index = position in `git log --first-parent e792cc40^..864a00d8` (161 commits).

| index | commit | date | test_heatNd_poincare | time |
|---|---|---|---|---|
| 1 | e792cc40 | 09/28 | pass, 44 checks | 43.4 s |
| 44 | 5308a063 | 10/20 | pass, 44 checks | 26.4 s |
| 87 | e56c7981 | 11/03 | pass, 44 checks | 28.9 s |
| 88 | e3a1f667 | 11/04 | error, affine 6.7e-01 | 24.4 s |
| 124 | f2f855d8 | 11/06 | error, affine 6.7e-01 | 22.1 s |
| 161 | 864a00d8 | 11/08 | error, affine 6.7e-01 | 23.5 s |

The diff of e3a1f667 names the equality writer, so it and its parent were measured first, then one
commit inside each half. That makes six runs, not a full binary bisection. As in a bisection,
behaviour between measured points is assumed monotone.

Between e56c7981 and e3a1f667 seven files change. On the test's path these are:

- `collect_eq_rows.m` – the new `drop_multiples`;
- `lpi_eq_cdopvar.m` – a new field of the row record only.

Byte-identical across the commit: `heatNd_poincare.m`, `heatNd_sdp.m`, the test, `sossolve.m`,
`lpi_soseq.m`, `impose_eq_rows.m`, `poscopvar.m`, `poscopvar_direct.m` and `copquadvar.m`.

## The failing family

The full error stack puts the failure in section (7), `test_heatNd_poincare.m:155`, not in
section (3).

- The family: 1-D, Dirichlet–Dirichlet, L2 form, tensor degree d = 0, face generators.
- It is the only family of the test whose basis misses target monomials. It has 5 equations with
  no decision variable (the 09/28 record: "uncovered rows 5").
- At e3a1f667 all 17 families of section (3) keep b affine (error ≤ 1.4e-16) – measured.

Measured with `heatNd_bisect_probe7.m` (the family rebuilt as `build_family` does, row record on):

- **e56c7981.** m = 10 at all three λ.
  - The five equations without a variable are rows 5, 6, 7, 9 and 10: cell 2, positions 6, 7, 8,
    13 and 14.
  - Their right-hand sides are b_6 = 1 + λ/3, b_7 = −λ/2, b_8 = λ/6, b_13 = −λ/6, b_14 = λ/6. These
    are fitted to the three measured values and exact to 7 digits.
  - b is affine to 5.9e-17.
- **e3a1f667.** m = 6 at all three λ, and one equation without a variable, row 5.
  - Written at λ = 0.5 μ = 4.935: position 6 (b = +2.645).
  - Written at 0.75 μ = 7.402 and 0.9 μ = 8.883: position 7 (b = −3.701, −4.441).
  - The row record lists the other four as dropped multiples of the written one, with factors that
    change with λ. For position 6 against position 7: −0.9369 at 0.75 μ, −0.8918 at 0.9 μ.
- **The switch.** |b_6| = |b_7| at λ = 6 (0.608 μ), between the first two build points.
  - Extrapolating row 5 from the first two points gives −7.509 at 0.9 μ; the written value is
    −4.441.
  - All other rows are affine, so the error is 3.067/‖b‖ = 3.067/4.5525 = 0.674 – the test's
    6.7e-01.

**Why At is unchanged.**

- The written equation of the group has no decision variable, so it is an all-zero column of At.
- The two members that win, positions 6 and 7, are adjacent among the written equations. So the
  zero column sits at index 5 whichever wins, and At is bitwise identical (`isequal`, measured).
- Had a written equation lain between the two winners, the test would have stopped one check
  earlier, with "At or K depends on lam".

## File and line

HEAD 8c416566; the code lines are identical since e3a1f667 (comments stripped and compared).

- `lpi_programming_sopvar/misc/private/collect_eq_rows.m`
  - :327 – row 1 of `Cs{k}` is the constant term A_k(c); `lpi_soseq.m:128` writes it as b.
  - :339 – `drop_multiples` runs over these columns, so b takes part in the grouping.
  - :452–457 – the threshold and the scale s = sign(first compared entry) × column maximum, both
    over all rows including row 1. An equation 0 = b_k has one compared entry (row 1), so s = b_k.
  - :473 – the key (count, first row, hashes). Every 0 = b_k has count 1, first row 1 and scaled
    entry ±1, so all of them share one key.
  - :479 – the representative is the largest |s|, i.e. the largest |b_k(λ)|. :519 drops the rest.
- `testing_codes/private/PIETOOLS_demos/sopvar_demos/heatNd_poincare.m`
  - :61–66 (header) – "Only b depends on lam". Every later λ uses b(λ) = b1 + (λ − λ1) db on the
    same At.
  - :280 – the build points 0.5, 0.75 and 0.9 μ.
  - :289–290 – the At and K identity check.
  - :295 – the affinity error.
  - :320, :323 – the interpolated b, and the deficiency flag computed from it.

## The guard is not a false alarm

Inferred, by arithmetic on the measured values:

- The interpolated row-5 right-hand side vanishes at λ ≈ 5.963 (0.604 μ).
- There every member of the group is nonzero (b_6 = 2.988, b_7 = −2.982).
- `sdp_at` would therefore report that program as not deficient.

## Masked failures in section (3)

Measured with a scratch copy of the test that catches only section (7). At e3a1f667 and at
864a00d8 the same nine checks fail:

| family | m written | formula / stock | dropped |
|---|---|---|---|
| N=2 L2 d=1 none | 481 | 512 | 31 |
| N=2 L2 d=1 product | 1057 | 1152 | 95 |
| N=2 L2 d=1 faces | 729 | 768 | 39 |
| N=2 L2 d=2 none | 2965 | 3042 | 77 |
| N=2 L2 d=2 product | 5029 | 5202 | 173 |
| N=2 L2 d=2 faces | 3889 | 3978 | 89 |
| N=2 energy sel=1 d=1 faces | 824 | 844 | 20 |
| N=2 energy sel=2 d=1 faces | 824 | 844 | 20 |
| (7) 1-D DD L2 d=0 | error | – | – |

- The dropped 2-D equations have decision variables. For those families b stays affine and At
  identical (measured).
- So these eight are expectation failures: the formulas count every nonzero equation, and the
  writer now omits verified multiples.
- The 1-D families drop nothing.
- Sections (1), (2), (4), (5), (6), (8) and (9) pass.
- The e3a1f667 commit message lists `test_heatNd_lpi` and `test_lpi_eq_sop` as run, not
  `test_heatNd_poincare`.

## Options

As written before the fix. The maintainer chose the library option and the section (3) counts
("Yes, do both"); see "Fix" above. The test-side option for (7) was not taken.

- **Library, `drop_multiples`.** Leave out of the grouping every column whose only compared entry
  is row 1 (count 1, first row 1: an equation 0 = b_k).
  - Effect (inferred): the written equations without a variable return to those of e56c7981, so (7)
    would again see 5 uncovered rows.
  - Cost (inferred, not measured): a mask on two arrays that already exist, O(number of equations).
  - The 3-D census in the e3a1f667 commit message needs re-measuring before adoption.
- **Test, `heatNd_poincare`.** Interpolate only equations with a variable, and compute deficiency at
  each λ from the program itself.
  - The library stays as it is.
  - The family then no longer reproduces the written program for deficient bases.
- **Either way, section (3) needs its counts changed.** Compare m plus the number of dropped
  equations (row record `prog.sopeq`, field `entries.dup`) with the formulas, or restate the
  formulas as written rows.

## Files

- In this folder:
  - `heatNd_bisect_check.m` – the runner: `heatNd_bisect_check(<export>/PIETOOLS, tag, logfile)`.
  - `heatNd_bisect_probe7.m` – the family probe: `heatNd_bisect_probe7(<export>/PIETOOLS, outfile)`.
- To recreate an export: `git archive <commit> PIETOOLS | tar -x -C <dir>`.
- Scratch logs and the exports were session-local; the exports are deleted.
