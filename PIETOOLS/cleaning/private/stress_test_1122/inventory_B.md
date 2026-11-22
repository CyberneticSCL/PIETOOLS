# Stress test of 11/22/2026 – inventory B (demos, example libraries, PIESIM, GUI files, root scripts, SOSTOOLS)

Queue: `inventory_B.tsv` (365 rows, 12 columns; 355 runnable, 10 not). Built read-only against HEAD 1bdb2114 of
branch ndopvar; nothing in the queue was run except six probes of under 30 s each that checked the call forms
(libcheck with a one-row case file, readtable of that file, a DDE block, a GUI file, PIESIM in 1-D and 2-D, the
library through the input shim). Generator and data (gen.pl, pde_rows.txt, rows_*.txt, blocklines.pl):
`C:\Users\mpeet\AppData\Local\Temp\claude\363c499c6332a503\874c85cb-c18d-4876-bea3-39cdd053e21f\scratchpad\inv_b\`.

## Runner contract (what the commands assume)

- One `matlab -batch` process per item; kill the process tree at `cap_s`. `est_s` includes about 15 s of start-up
  and path setup per process (night runner of 09/29: 17 s wall for items of 1–3 s).
- Path: `restoredefaultpath` first (the saved pathdef points at the main checkout), then the worktree tree as
  `lb_setup.m` sets it (heatNd_path: one addpath(genpath), nested worktrees stripped, single resolution asserted),
  then MOSEK 11 (`C:\Program Files\Mosek\11.0\toolbox\r2019b`) and SeDuMi 1.3 (`C:\YALMIP\yalmip_install\tbxmanager\toolboxes\sedumi\1.3\pcwin64`).
  SeDuMi is needed by the SOSTOOLS demos and by PIETOOLS_PDE.m (solver 'sedumi').
- `cd(fullfile(<worktree>,'PIETOOLS',folder))`; `set(0,'DefaultFigureVisible','off')`; `rng(1)`; evaluate `cmd` in
  the base workspace (`evalin('base',cmd)`: scripts call clear, libcheck clears the base workspace); `close all force`.
- Prompts. `input()` errors under -batch (MATLAB:services:MissingRequiredCapability, probe 11/22). Items that pass
  through an `input()` prompt put `cleaning/private/stress_test_1122/shim` first on the path (absolute path of the
  main checkout; the folder is untracked, so a worktree does not contain it) and set the answers with
  `setenv('ST_INPUT',...)`, '|' between answers. Used by PDE_16, PDELIB_EXTRACT, PDELIB_BATCH (n), ROOT_PIETOOLS_PDE (y)
  and PIESOS_LIB_* (y, then Enter: the settings of the example).
- Output markers: libcheck prints `##### case <id> <executive>: <status> <value>`; the other items print `## ST ...`
  lines or the executive summaries (`[stability_sop P] optimal ...`).

## Routes

- PDE library (84 runnable): the acceptance harness libcheck.m (testing_codes/private/baselines/2026-11-18_studies/private,
  tracked) with a one-row case file written to tempdir: the example file as the library calls it, its st, its flags and
  forms, PIETOOLS_auto_execute_sop. The library function itself (switch, both output formats, prompt) runs in two loop
  items through the shim; the batch-format item converts each batch PDE and solves light stability.
- DDE and NDS/DDF libraries (39 runnable): the library is copied to tempdir with one block uncommented (its code lines
  and its st and form lines, line numbers of HEAD; each command asserts the line count 914 / 579) and put first on the
  path, then PIETOOLS_DDE.m or PIETOOLS_DDF.m runs; PIESIM is stubbed for the DDE blocks, as blockcheck.m does
  (PIETOOLS_DDE.m simulates with uinput.ic = [5,5,5] before the LPIs). The active blocks (DDE26, DDF03) run in
  ROOT_PIETOOLS_DDE (with PIESIM) and ROOT_PIETOOLS_DDF. blockcheck.m is not used: it extends a block only over lines
  tagged 'MMP, 11/19/2026', and the comment pass of 11/21 retagged every st line (no such tag is left in either
  library), so at HEAD it would run the blocks without their st and forms (read from the code, not run; a defect
  of the harness, which inventory A owns).
- GUI save files (38): parsePDEParams on a hidden uitable, convert_PIETOOLS_PDE_terms (legacy terms converter), light
  stability and, where w and z exist, light H-infinity gain. Records of 11/08 predate the tier change of 11/17.
- PIESIM: examples_pde_library_PIESIM_1D/2D at the solver_PIESIM defaults (N 16, tf 0.1, dt 1e-3, BDF2,
  ifexact true), plot and ploteig 'no'.

## Counts and estimated time (runnable items)

| Group | Rows | Runnable | Est. (s) |
|---|---|---|---|
| DEMO (DEMO1–15; heatNd_apply) | 16 | 15 | 570 |
| SNIPPET (snippets_from_manual) | 9 | 9 | 1,005 |
| PDE_LIB (cases 1–72, variants, files X1–X4) | 87 | 84 | 7,315 |
| PDE_LIB_ROUTE (library function) | 4 | 2 | 330 |
| CLI (command_line_input_Examples) | 33 | 33 | 660 |
| GUI (*_GUI.mat, Demos/*.mat, app) | 39 | 38 | 950 |
| DDE_LIB | 26 | 26 | 550 |
| NDSDDF_LIB | 15 | 13 | 260 |
| FI_LIB (17 inequalities) | 17 | 17 | 340 |
| PIESIM_1D (1–41) | 41 | 41 | 820 |
| PIESIM_2D (1–32) | 32 | 32 | 1,120 |
| ROOT (PIETOOLS_PDE/DDE/DDF/FI, solver_PIESIM, path script, installer) | 7 | 6 | 220 |
| SOSTOOLS (22 demos, addsostools) | 23 | 23 | 3,175 |
| PIESOS (3 demos, 13 library examples) | 16 | 16 | 1,925 |
| **Total** | **365** | **355** | **19,240 (5.3 h)** |

Sum of caps: 130,965 s (36 h), the worst case. Estimates marked "inferred" in `source` have no recorded run:
the PIESOS items, DEMO14/15, the snippets Ch2/Ch6/Ch11, sosdemo11s/12s, PDE_70, PIESIM at N 16.

## Longest items and memory

| Item | est / cap (s) | Peak memory (GB) | Note |
|---|---|---|---|
| PDE_70 (2D_Wave_Acoustic_Boundary_Beale) | 2,715 / 8,145 | 30 (model) | light Q, 2,014,033 vars, m 30,624, never solved; run alone, last |
| SOS_sosdemo12s | 1,800 / 5,400 | 4 | timed out at 3,600 s on 09/29 (8-core low-priority mode) |
| PDE_54 (3D_Advection_Diffusion_Eq) | 1,150 / 3,450 | 15 (inferred) | 1,835,012 vars, m 26,429; 664 s (11/19), 1,134 s (11/21) |
| SOS_sosdemo11s | 600 / 1,800 | 2 | never timed |
| PDE_50 (PIETOOLS_PDDE_Ex) | 595 / 1,785 | 10 (inferred) | H2 503,633 vars, dual 1.66e6 vars; 577 s |
| PDE_36 (2D_Telegraph_Eq) | 520 / 1,560 | 18 (model) | light Q, 1,775,200 vars, m 25,947; 501 s (11/18) |
| PDE_39 (2D_NS_Antonelli) | 280 / 840 | 25 (model) | light Q, 648,900 vars, m 30,721; 260 s |

Model: MOSEK peak about 27.1 GiB x (m/32,256)^1.89, which overstates 3-D programs about twice. The five items of
10 GB or more (PDE_70, 39, 36, 54, 50) need the machine to themselves: no second MATLAB process, a free-memory guard.

## Not runnable (10)

| Item | Reason |
|---|---|
| PDE_55 (3D_Coupled_Heat_Eq) | file st: 7,340,048 vars, m 105,717; Sol job cancelled at 6 h, 191 GiB; needs a Sol high-memory node |
| PDE_56 (2D_Delayed_Heat_Eq_Fridman) | file st: 2,524,300 vars, m 51,558, above the 64 GB MOSEK limit (near m 45,000) |
| PDE_69 (Wave_Heat_Transmission_2D_Batty) | file st: 4,453,560 vars, m 60,678, above the 64 GB limit |
| PDELIB_GUI_OPTION | library option 'gui' opens the GUI window (the save files run headless in group GUI) |
| PDELIB_100_101 | cases 100, 101 cannot be selected (indices up to 72); 101 reads undefined A, B, C |
| GUI_APP | PIETOOLS_PDE_GUI.mlapp is interactive (its conversion calls run in the GUI items) |
| ROOT_pietools_update | downloads from GitHub, installs tbxmanager and SeDuMi, edits startup.m, saves the path, prompts |
| SOPVAR_heatNd_apply | helper function with required inputs; the heatNd benchmark scripts are in testing_codes/private |
| NDS06 | block disabled in the library (two comment levels) |
| NDS_FRAG | fragment (An never defined) |

## Expected failures recorded earlier (40 runnable items carry a note)

- Not certified as the files stand: PDE_36, 39, 41 (unknown), PDE_44, X3 (infeasible), PDE_50 dual cco
  (infeasible), PDE_70 (infeasible at every setting tried); PDE_15 prints optimal but fails the gate of 11/20.
- PDE_16: error by design. DDE_RAND: expected error, the pvar s is undeclared (inferred). DDE08: the block now sets
  form Q after P, so Q/Qd run against a P record. PIESIM1D_32: error at PIESIM_int_bdf:289 (11/20).
- Status unknown since 09/29: CLI_stability_8, CLI_example_2D, SNIP_Ch2 (errors on 09/29).
- Seven GUI save files differ from their library files and four by rounding (data, 11/05); DEMO13 prints
  "certifies no bound" (11/18).

## Scope notes

Added beyond the brief, since no other inventory covers them: PIESOS demos and library, the 33 command-line example
scripts, solver_PIESIM.m, pietools_path_update.m, addsostools.m, the example files without a library index
(Reaction_Diffusion_Robin, Wave_Eq_w_Boundary_Delay_2D, 2D_Plate_Eq, 2D_Wave_Eq_Damped; the last two only convert:
they set no executive flag by design). Left out as not user-facing: PIETOOLS_demos/cuadmm and lowrank_2d_stability
(brief), orphaned_file/private (10 GUI files, old snippets, ndopvar tests, stability_16*.m), PIETOOLS_demos/html,
Contents.m files, the editor autosave examples_DDE_library_PIETOOLS.asv, and testing_codes (inventory A, which holds
test_examples_FI_library.m).

Overlap with inventory A: its ST_libcheck runs libcheck_cases.csv in one process (71 of the 84 PDE_LIB cases, no
cap per case) and its ST_BLOCKCHECK runs blockcheck.m, which at HEAD drops the block st (see Routes). One route per
example is enough; the per-case and per-block items here isolate failures and caps.
