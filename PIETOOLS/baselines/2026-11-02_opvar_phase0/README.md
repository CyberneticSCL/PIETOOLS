# Phase 0 of the opvar removal: baseline fixtures, 11/02/2026 and 11/04/2026

Reference results recorded while the classes `opvar`, `dopvar`, `opvar2d` and `dopvar2d` still exist,
so that each later phase of the removal can be checked against data the removed code produced.
The plan is `opvar_removal_plan.md` (Phase 0, §3.6, §4).

## Tree

| Passes | Commit | Export | Recorded |
|---|---|---|---|
| 1, 2 | `eb5d0600` (merge of 10/30) | `C:\Users\mpeet\wt\p0_eb5d0600` | 11/02 |
| 3, 4 | `c81c3d04` (ndopvar HEAD of 11/04) | `C:\Users\mpeet\wt\p0_c81c3d04` | 11/04 |

**Passes 3 and 4 are the reference for the later phases.** Passes 1 and 2 stay, so that
`compare_pass1_pass3.txt` shows what the commits between the two trees changed.

- Every run uses a frozen export of the commit, made with `git archive <commit>:PIETOOLS`
  (`P0_COMMIT.txt` names the commit). The environment variable `P0_ROOT` names the export; the
  default is `p0_c81c3d04`. Nothing of the main checkout is on the path (`private/p0_path.m`; it
  asserts that every name the fixtures use resolves to exactly one file inside the export).
- Each export is deleted after its runs. To record another pass on one, recreate it from the
  repository root (`GitHub/PIETOOLS`), here for `c81c3d04`:

  ```
  mkdir -p /c/Users/mpeet/wt/p0_c81c3d04
  git archive c81c3d04:PIETOOLS | tar -x -C /c/Users/mpeet/wt/p0_c81c3d04
  git rev-parse c81c3d04 > /c/Users/mpeet/wt/p0_c81c3d04/P0_COMMIT.txt
  ```
- Between the two trees (ten commits of 11/01–11/04): the equality-row dedupe of the LPI layer
  (`e3a1f667`: `collect_eq_rows`, `lpi_eq_cdopvar`, `lpi_eq_sdopvar`), the optional cuADMM solver and its
  classifier status (`c1374587`, `d4ca0722`; never the default), and Declan's converter work: the N-D
  converter (`b4256de6`: `convert_PIETOOLS_PDE`, `initialize_PIETOOLS_PIE`, `expand_tderivatives`,
  `@pde_struct/convert`, `parse_copvar_spaces`, `mat2copvar_grid`, `DEMO10_PBC`; the BVP solver of
  `e53d8801`), new `sopvar` methods and `@sopvar/mtimes` (`bae0f3f6`), PIE/NPIE converters and
  polynomial addition across output spaces (`0d537053`, `e56c7981`), and the 3-D library examples 53–55
  (`6e65e806`). The library change only appends cases 53–55, so catalog rows 1–52 build the same
  examples; the 10/28 catalog, and with it F2 and F3, does not contain 53–55.
- The `eb5d0600` export does **not** contain the uncommitted edits of the main checkout at the time
  (`SOSTOOLS400/sossolve.m`, `executives_sopvar/utility_functions/exec_info_sop.m`,
  `lpi_solve_loop_sop.m`, `lpi_programming_sopvar/misc/lpi_classify_sop.m`, `lpi_eq_cdopvar.m`,
  `lpi_eq_sdopvar.m`, `misc/private/collect_eq_rows.m`, …), among them the equality-row change that
  lowers the row count m. F3 and F4 depend on those files: re-record them (a new pass) on any tree
  whose LPI layer differs. F1 and F2 do not call them.
- MATLAB R2025b, MOSEK 11.0.30.

## Contents

| Fixture | What | Files |
|---|---|---|
| F1 | PIESIM on the 41 1-D and 32 2-D library examples (`examples_pde_library_PIESIM_1D/2D`), N 16 (1-D) or 8 (2-D), tf 0.1, dt 1e-3, BDF 2: the solution and grid PIESIM returns, and every discretized operator PIESIM integrates (`Dop`), captured by a shadow of `PIESIM_time_integrate` that forwards to the original; for 1-D examples with an exact solution, the error at the final time as `PIESIM_plot_solution` forms it. | `F1/pass<k>/<d>D_ex<NN>.mat`, `index.csv` |
| F2 | Every library PIE: the 41 + 32 PIESIM examples and one row per example of the 10/28 catalog (PDE, DDE, NDS/DDF libraries). Per nonempty operator (T, Tw, Tu, A, Bw, Bu, Cz, Dzw, Dzu, Cy, Dyw, Dyu): the operator as a `copvar`, 3 random polynomial test functions per input component, and the output at 5 random points per output component, by `apply_sopvar` block by block. For PIEs that hold opvar operators (the legacy input formats), also `apply_opvar` on the opvar at the same points, and the opvar fields as a plain struct readable without the class. | `F2/pass<k>/<label>.mat`, `index.csv` |
| F3 | The container (`_sop`) executive rows of the 10/28 catalog through the copied 10/28 harness (`private/lb_run.m`, path `sop`), tier light (heavy not recorded, see Results); each row's note ends with the classifier status, the SDP shape (m, nx, ndv, Ks) and the witness. PDE 36 and 38 are skipped (1,653 s and 1,638 s at light on 10/28). PDE 50 and 50/H2_norm_dual are skipped from 11/02 06:30 on: an error on 10/28, on the export one light solve of the primal took 2,421 s; the pass 1 light primal row, recorded before the skip, stays. | `F3/pass<k>/sop_results.csv` |
| F4 | DEMO1–DEMO11, each in its own MATLAB process: the script's variables after it ends (LPI programs dropped) and its time. | `F4/pass<k>/<DEMO>.mat`, `index.csv` |

Passes 1 and 2 are on `eb5d0600`, passes 3 and 4 on `c81c3d04`; pass 3 is the reference for the
later phases. The second pass on each tree repeats F1, F2 and F3 (light); F4 has one pass per tree.
`compare_pass<a>_pass<b>.txt` gives the differences of the stored outputs (`private/p0_compare.m`),
`compare_F2_pass<a>_pass<b>.txt` those of the F2 operators independent of block partition and state
order (`private/p0_f2cmp.m`).

## Running

```
bash run_phase0.sh now           # eb5d0600: F1, F2 pass 1
bash run_phase0.sh after0430     # eb5d0600: F3, F4, pass 2, comparison; waits until 04:30
bash run_phase0.sh rest          # eb5d0600: F3 light rows not yet recorded, both passes; comparison
bash run_phase0.sh c81now        # c81c3d04: F1, F2 pass 3; comparison with pass 1
bash run_phase0.sh c81after0430  # c81c3d04: F3 light passes 3 and 4, F4 pass 3; comparisons
```

F1 and F2 pass 4 were run by hand on 11/04 (01:54–02:00, `p0_f1(4)`, `p0_f2(4)` with `P0_ROOT` set;
logged in `runs_phase0.log`), and `p0_f2cmp(1,3)` likewise.

Each step is one MATLAB process (`matlab -sd private -batch ...`) with a time cap; the log is
`runs_phase0.log`, the MATLAB output `logs/<step>.txt`. A file `STOP` in this folder skips the steps
not yet started.

## Use in later phases

Record a new pass of the same fixture on the changed tree and compare it with pass 3: export the
tree, set `P0_ROOT` to it, run `p0_f1(k)`, `p0_f2(k)`, `p0_f3(k,{...})` or `p0_f4one(name,k)`, then
`p0_compare(3,k)` and `p0_f2cmp(3,k)`. Two passes on one tree bound the run-to-run variation; a
difference above it belongs to the change. After the legacy input formats are retired, the opvar
outputs (`Yop`) of F2 are the reference for the same examples posed in `pde_var` form, and
`p0_compare` reports `Yop` as different by design.

`p0_compare` compares the stored F2 outputs, which `p0_apply` draws per block: a converter that splits
a container into other blocks, or orders the state components differently, changes them although
the operators are equal (measured 11/04, below). For F2, `p0_f2cmp` is the comparison: it applies
both passes to the same test functions drawn per scalar component, and for a PIE that still differs,
`p0_f2perm` searches one permutation of the state components for all of its operators.

No fixture file stores an `opvar`, `dopvar`, `opvar2d` or `dopvar2d` object (measured 11/02 by loading
every file and walking its structs, cells and properties, `p0_classes`): F1 holds `pde_struct`, `polynomial` and
`sym` objects, F2 `copvar` and `polynomial`, F4 also `cdopvar`, `dpvar` and `pie_struct`; the opvar
data of F2 is a plain struct. The fixtures therefore load after the four classes are removed, provided
these classes are on the path. With MATLAB's factory path, a load of `F4/pass1/DEMO6_*.mat` returned
its `pie_struct`, `copvar` and `polynomial` variables as 0×0 doubles and printed no warning. The saved
MATLAB path of this workstation puts the main checkout's `opvar/` on the path but not `sopvar/`.

## Files

- `private/`: the harness. `p0_*.m` are new; `lb_*.m`, `catalog.mat` and `blocks/` are copies of
  `baselines/2026-10-28_example_library/private`, with three changes marked "Phase 0, 11/02/2026":
  `lb_setup` (the export's path), `lb_build_system` (PIE size of a copvar `T`), `lb_exec` and
  `lb_run` (the executive's info output in the note).
- `private/capture/PIESIM_time_integrate.m`: the capture shadow (on the path only during F1).
- `private/p0_classes.m`: `p0_classes('F4/pass1')` lists the classes of the objects stored under a
  fixture folder, to check that the fixtures still load on a tree without the four classes.
- `private/p0_f2cmp.m`, `private/p0_f2perm.m`: the F2 comparison by components and up to a state
  permutation.

## Results

**Pass 1, F1 and F2 (11/02, 01:56–02:08).**

- F1: 72 of 73 examples recorded (7.2 min). 1-D example 32 fails inside PIESIM with these options
  (`PIESIM_int_bdf:237`, incompatible array sizes); the error is the fixture. The 38 legacy-format 1-D
  examples take the batch branch of `convert_PIETOOLS_PDE` (16) or the legacy-terms branch (22), and
  PIESIM receives their PIE with opvar operators; the other 35 examples arrive as copvar PIEs and pass
  through `pie2opvar`. 6.2 MB.
- F2: 172 PIEs recorded (4.4 min): 73 PIESIM examples and 99 catalog examples; 38 hold opvar operators,
  134 copvar; 47 have two spatial variables; 1,134 blocks applied. On the 38 opvar PIEs, `apply_opvar`
  and `apply_sopvar` of `op2copvar_sop` agree to a relative difference of at most 2.6e-14 (median
  1.2e-16). 1.3 MB.
- Smoke test of F3 and F4 before queuing (outputs deleted): PDE 1 stability, light, certified, m 90,
  Ks [14 10 7], maxre −8.852; PDE 20 H∞ gain, light, γ 0.50081445 with witness gain 0.5; DEMO2,
  operator_norm 0.6368044657 in 1 s.

**Pass 1, F3 and F4; pass 2 and the comparison (11/02, 04:30–06:27 and 06:30–06:39).**

- F3, light: 135 rows (69 PDE, 41 DDE, 25 NDS/DDF). Status: certified 89, converged 20, unknown 15,
  infeasible 5, hi_pass 2, lo_fail 2, error 2 (DDE19 and DDF03, `PIETOOLS_Hinf_control_sop`,
  `parse_copvar_spaces:158`; both also error on 10/28). Against the 10/28 light rows: 132 rows
  have the same status, of which 128 have the same result (relative difference at most 1e-6) and 4
  have NaN in both; the other 3 errored on 10/28 and solve on the export: PDE 42 certified,
  PDE 43 infeasible, PDE 50 unknown.
- PDE 50 (`PIETOOLS_PDDE_Ex`, H2 norm, 2-D) errored in 4.7 s on 10/28. On the export the light primal
  solves in 2,421 s: ndv 5,080,609, m 40,881, 13 PSD blocks of order up to 1,261, status unknown
  (numerr 2), result 263.6 against the documented 0.8502 (custom settings with psatz). The 4,500-s
  cap of pass 1 stopped its dual after 23 min, and the 1,800-s cap of pass 2 stopped the primal after
  19 min. Both rows are skipped from 06:30 on (`p0_f3`); the other light rows ran in 8.5 min
  (`run_phase0.sh rest`). PDE 50 took 82 of the 126 min of runs after 04:30. The 75-min estimate used
  the 10/28 times, in which PDE 50 errored in 5 s.
- F3, heavy: not recorded; the pass 1 cap expired inside the light tier. Estimate: 29 min per pass,
  from the 10/28 heavy times of the rows other than PDE 36, 38 and 50 (1,309 s) times 1.31, the ratio
  of the light times of the same 66 PDE rows (707.5 s on the export, 541.5 s on 10/28). The heavy
  settings of 2-D rows changed in `e78ed85b`, so the rows the projection skips at heavy can differ
  from 10/28.
- F4: DEMO1–DEMO11 ran to the end, 141 s in all (`F4/pass1/index.csv`): DEMO2 operator norm
  0.6368044657; DEMO3 Poincaré constant 0.3876803676; DEMO4 lam 9.8828125 after 8 iterations
  (lam_min 9.84375, lam_max 9.921875); DEMO5 γ 1.000093476; DEMO6 γ 1.414977876; DEMO7 γ 1.000126864
  (control) and 2.758965915e-4 (observer); DEMO8 H2 norm 0.101574958; DEMO9 γ 0.09397896859; DEMO10
  lam 8.869604401; DEMO11 lam 20.7392088. DEMO1 has no scalar result other than its parameters
  (b 0.01, c 1); its file holds 32 variables, among them its PIE (`pie`) and PIESIM solution
  (`solution`). DEMO11 warns at
  `getController_sop:148` that its gain K depends on 2 spatial variables and has no opvar form
  (`Kop = []`); DEMO7 warns in PIESIM that it has more initial conditions than PIE states.
- Pass 2 against pass 1 (`compare_pass1_pass2.txt`): F1, 73 examples, and F2, 172 PIEs, identical
  (largest relative difference 0 in the captured operators, the solutions and the operator outputs);
  F3, 134 light rows in both passes, none with a different status, result or SDP shape; PDE 50 light
  in pass 1 only. F4 has one pass.
- Not recorded: F3 heavy (both passes), the PDE 50 dual, and PDE 50 in pass 2. Each needs a night
  window (heavy about 29 min per pass, PDE 50 about 40 min per solve) and the export recreated; for
  PDE 50, remove its two rows from the skip list in `p0_f3`.

**`c81c3d04`, F1 and F2, pass 3 against pass 1 (11/04, 01:39–01:45; 5.2 min).**

- F1 (`compare_pass1_pass3.txt`): 48 of 73 examples identical in `Dop` and the solution, among them
  39 of the 41 1-D examples; 1-D example 32 errors in both passes, 1-D example 41 differs by 2.4e-10
  in `Dop`. 21 2-D examples differ by at most 1.1e-9 in `Dop` and 1.3e-9 in the solution. 2-D example
  08: `Dop` 2.0e-8 (the inverse of `Tcheb` has entries of 1.6e9), solution 0.34, the solution reaching
  |x| = 4.7e38 by t = 0.1, so that the change of `Dop` grows with it; 2-D example 16: `Dop` 8.8e-8,
  solution 2.5e-4 (|x| ≤ 31). Pass 1 and pass 2 on one tree were identical, so these differences come
  from the commits between the trees.
- F2, stored outputs (`compare_pass1_pass3.txt`): 102 of 172 PIEs identical, 40 differ by at most
  1.3e-14, 7 by 1.5 to 15.3, and 23 have outputs of other sizes. The converter of `b4256de6` builds
  other blocks (58 PIEs changed their number of blocks, e.g. a 4-component state in one block now in
  four), and `p0_apply` draws its test functions per block, so these numbers do not compare operators.
- F2 by components (`compare_F2_pass1_pass3.txt`, `p0_f2cmp`): 641 of 668 operators agree to at most
  1e-12. The other 27, in 10 PIEs (catalog PDE 15, 16, 19, 32, 39, 41, 44; PIESIM 2-D 19, 20, 28),
  agree to at most 3.5e-12 after one permutation of the state components per PIE, applied to all of
  its operators (for PIESIM 2-D 28 the state order is unchanged and the disturbance components are
  permuted). Apart from block partition and component order, no F2 operator changed.
- Pass 4 against pass 3 (11/04, 01:54–02:00): F1, 73 examples, and F2, 172 PIEs, identical (largest
  relative difference 0).

**`c81c3d04`, F3 and F4 (11/04, 02:04–02:35; started at once at the maintainer's request, not at
04:30).** Another session's MATLAB process (started 02:00, 4.2 GB) ran at the same time, so these
times are not comparable with 11/02; F3 pass 3 took 13.5 min, pass 4 12.2 min, F4 4.4 min.

- F3, pass 4 against pass 3 (`compare_pass3_pass4.txt`): 134 rows, none with a different status,
  result or SDP shape.
- F3, pass 3 against pass 1 (`compare_pass1_pass3.txt`): 134 rows in both (PDE 50 only in pass 1).
  - m: lower on 116 of the 132 rows with an SDP in both passes, higher on none; 154,978 → 148,863 rows
    over those 116 (−3.9%), at most −29.0% (PDE 48, 5,298 → 3,764).
  - nx, ndv or Ks: changed on PDE 15, 16, 18, 19, 44 (e.g. PDE 15, ndv 21,888 → 19,188, Ks
    [112 84 40] → [104 78 40]); their status and result did not change. These are PIEs the converter
    of `b4256de6` splits into more blocks; the positive-operator construction then builds smaller
    Gram bases.
  - Status, 5 rows: PDE 3 unknown → certified; NDS09/stability_dual unknown → certified;
    NDS11/stability_dual certified → unknown; DDE19 and DDF03 (`PIETOOLS_Hinf_control_sop`) error →
    unknown and inaccurate, the error at `parse_copvar_spaces:158` gone (that file changed in
    `b4256de6`).
  - Result, 17 more rows with the same status: certified H∞ bounds by 1.8e-6 to 1.1e-3 relative
    (PDE 23, 24, 26; DDE11, DDE12, DDE12 dual, DDE13, DDE18, DDE26 dual; largest DDE18 0.3791484 →
    0.37801951, PDE 23 23.709405 → 23.718499), bisection end points DDE01 2.0160938 → 2.0180859 and
    NDS08 0.027307861 → 0.027523096, and the values of 6 rows with status unknown (PDE 30, 45, 46, 47,
    48, DDE26).
  - Every one of these 22 rows has a PIE whose operators agree with pass 1 to at most 1e-12 by
    components (`compare_F2_pass1_pass3.txt`; for 15 of them the stored F2 outputs are bit-identical).
    None of the 10 reordered PIEs changed status or result. The changes therefore come from the LPI
    layer (the equality-row dedupe of `e3a1f667`) and, for DDE19 and DDF03, from `parse_copvar_spaces`.
- F4, pass 3 against pass 1: all 11 DEMOs ran to the end; DEMO1, 7, 8 and 9 have no numeric
  variable that differs by more than 1e-6. Reported values: DEMO2, 8, 10, 11 equal;
  DEMO3 Poincaré constant, DEMO6 and DEMO9 γ within 7e-8 relative; DEMO5 γ 1.000093476 → 1.000094369
  (8.9e-7); DEMO7 observer γ 2.758965915e-4 → 2.759920741e-4 (3.5e-4 relative, 9.5e-8 absolute);
  DEMO4 lam, lam_min, lam_max and the 8 iterations equal, its last step now reported primal
  infeasible (`is_pinf` 0 → 1). Other numeric variables: the SDP dual of DEMO3 (3.3e-3), the
  closed-loop trajectories of DEMO6 (2.1e-6) and the estimator trajectories of DEMO5, which reach
  9.7e5 and differ by O(1) relative, as the computed estimator gain differs; 70 of 184 object
  variables (PDE, PIE, decision variables) differ.
