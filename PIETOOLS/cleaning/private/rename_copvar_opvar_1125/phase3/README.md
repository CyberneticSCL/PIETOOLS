# Phase 3 – verification, 11/27/2026

Identity before and after the rename, on exports of 4ac5c66d (parent of the rename) and abfb0e83
(rename and Phase 1b), outside OneDrive: `C:/wt1127/old`, `C:/wt1127/new` (`git archive`
without `cleaning/`, `orphaned_file/`, `.claude/`). A rename changes no arithmetic, so every
number must match.

## Method

- Runner `d3run.m`: `strun1122.m` (the 24-row runner of 11/22–11/24) adapted. Per row the path
  holds one export only, checked before and after the row (nothing of the checkout, no
  `cleaning/` or `orphaned_file/` code; the core routines and the container class – `copvar`
  before, `opvar` after – resolve inside the export); `rng(1)`; figures off; the capture folder
  first on the path.
- Capture (`capture/`): shadows of `Sedumi2Mosek` (A = At', b, c, K: the SDP that `sossolve`
  hands to MOSEK), `MosekSol2SedumiSol` (x, Y; `sossolve` takes y = Y(1:m)) and `sedumi`
  (inputs and outputs). Each forwards to the original through a handle bound with the capture
  folder off the path, and logs a SHA-256 digest of the exact contents of each array (class,
  size, sparse indices and values): equal digests mean equal bit for bit. Test instrument only;
  no toolbox file changes.
- Queue (`q_old.tsv`, `q_new.tsv`, the same row ids): T1 converter tests, manual snippets Ch3,
  Ch5, Ch7, `test_executives_sop`, 11 container LPI tests; T2 the class tests (`Test_copvar` /
  `Test_opvar` and the like), 7 `sopvar` tests, the heatNd PIE test, the FI library, a PIESOS demo
  and 5 PIESOS tests; D DEMO1–DEMO15, each base workspace saved class-free (`d3_ws.m`); P the
  phase-0 fixtures F1 (PIESIM on the 1-D and 2-D library examples) and F2 (every library PIE,
  each operator applied to random polynomials at random points); X the 3-D heat LPI
  (`e2e_heat3d_stability`).
- One MATLAB process per batch, which takes no new row after 17 min (`p3_drive.sh` resumes it):
  no run exceeds 20 min.
- Repeatability: the same SeDuMi row and MOSEK row run twice on the new tree give equal digests
  of A, b, c, K, x and y, so a difference across the trees would be a change, not noise.

Incident: the first run failed 8 rows in the harness – `d3_hash` passed an empty array to
`MessageDigest.update` for a sparse input with no nonzero (c = 0 of a feasibility SDP), and
MATLAB chose `update(byte)`. Fixed (class and size only for an empty sparse input), checked on
empty and nonempty inputs, and rerun from scratch; the failed run is discarded.

## Results

- Rows (`cmp_sdp_all.txt`, `results_old.tsv`, `results_new.tsv`): 55 rows ran on both trees
  (T1 18, T2 19, D 15, P 2, X 1), all OK on both. One run on the new tree ended MATLAB 2 s into
  `sv_test_sopvar_kron`: "Access violation", thread 'LXE compiler thread', no PIETOOLS frame
  (`matlab_crash_dump.109016-1`), the intermittent fault of the 11/22 stress test; rerun alone,
  the row passed (48.6 s; 46.8 s on the old tree). The D6 test (`test_call2025_sop`, new tree
  only) passed.
- SDPs: every row has the same sequence of SDP records on both trees – 110 SDPs equal bit for
  bit in A, b, c and K, 109 of them also in x and y (219 records): the converter and container
  LPI tests, the executives test, the manual snippet Ch7 (SeDuMi), the FI library (15), the
  PIESOS demo and tests, DEMO1–DEMO15 (DEMO4: 8 bisection steps), and the 3-D heat LPI (A of
  22,226 x 2,713,563 with 14,431,582 nonzeros, 1,100 and 1,090 s per row). The 3-D solve reads
  the MOSEK result itself (`heatNd_solve`), so its x and y are not captured; the test checks
  them and passed on both trees.
- DEMO workspaces (`cmp_ws.txt`): 15 files, no differing leaf except three wall-clock times in
  DEMO13 (`info_cl.hist.t_build`, `t_solve`, `witness.t`; 0.09 to 0.28 s), which the time-field
  filter did not name.
- F1 (`cmp_f12.txt`): PIESIM on 73 library examples (41 1-D, 32 2-D); the captured discretized
  operators and the solutions differ by 0 (largest relative difference).
- F2: 172 library PIEs, each operator applied to random polynomials at random points; outputs
  differ by 0.
- Cost (`classops.txt`, `t3_classops.m`): medians of six container operations in 1-D (R x
  L2^2[s]) and 2-D (R x L2[s1] x L2[s2] x L2[s1,s2]), two runs per tree, alternating: the trees
  overlap within the run-to-run spread (2-D T'\*P\*T: 228 and 243 ms before, 222 and 226 ms
  after; 1-D positive declaration: 34 and 44 ms before, 47 and 38 ms after), and the
  decision-variable counts are equal (1-D 1,035 and 97, 2-D 10,440 and 225). With identical SDP
  data, the cost in the decision variables and in the number of spatial variables is unchanged.

## The 2025 field assignments of the manual (Chapter 3)

The Chapter 3 snippet of v5.0 (`ch3_2025.m`), cell by cell on each tree (`ch3_record.m`):

| Section | Before (4ac5c66d) | After (abfb0e83) |
|---|---|---|
| 3.1.1 3-PI operators: `opvar A; A.I = [-1,1]`; `opvar B; B.var1 = s; B.var2 = s_dum; B.R.R0 = R0; ...` | Undefined function 'opvar' | runs: A an `opvar` 2 x 2, B an `opvar` 2 x 1 (containers; the 1-D fields are emulated) |
| 3.1.4 4-PI operator, dummy variable `tt` | Undefined function 'opvar' | error: "Dummy variables must have the default names, the spatial variables with '_dum' appended: s_dum. Set var1 first to use other names." |
| 3.2.1, 3.2.2 9-PI and general 2-D (`opvar2d`) | Undefined function 'opvar2d' | the same: `opvar2d` retired on 11/08 |

## Other checks

- SOSTOOLS400: no file differs between 4ac5c66d and abfb0e83. `sopvar/@sopvar`, `@sdopvar`: 12
  files, 19 lines, names and comments only (`dpvar_op_copvar` → `dpvar_op_opvar`).
