# Remove the `_sop` suffix – plan

Maintainer, 11/25/2026: "compile a list to functions with the cop/copvar suffix with a goal to remove
these suffixes to try and match the 2025 release naming conventions where possible"; then "Also
applies to subroutines". The container rename (`rename_copvar_opvar_1125/`) left the `_sop` suffix
for a later pass (its PLAN.md, Section 8). Maintainer, 11/28/2026, on that section: "Lets do this
now".

Status: decisions E1–E4 taken on 11/28 (Section 5). Phases 0 and 1 done: the rename committed as
fe6d3530 (`phase1/README.md`); Phase 1b 0d02ac08 (`phase1b/README.md`); Phase 2 182a2c49
(`phase2/README.md`); Phase 3 run 11/28: 56 rows, every SDP equal bit for bit, F1/F2 differences
0 (`phase3/README.md`); announcement drafted (`ANNOUNCEMENT.md`). Inventory at 91ea628e (=
origin/ndopvar). Table:
`rename_table.tsv`; scripts and raw outputs: `inventory/`.

## 1. Scope

- **Functions:** every toolbox function whose name ends in `_sop` (35) loses the suffix, taking the
  2025 name where one exists. The six tiers `settings_sopvar_<tier>` take the 2025 names
  `settings_PIETOOLS_<tier>`; `eq_opts_sopvar` (it serves `dopvar` and `sdopvar`) becomes
  `get_eq_opts` (2025: `get_eq_opts_2D`). 42 functions in all.
- **Folders:** `executives_sopvar`, `lpi_programming_sopvar`, `settings_sopvar` become
  `executives`, `lpi_programming`, `settings` (2025 names), with the test mirrors
  `testing_codes/private/executives_sopvar`, `.../lpi_programming_sopvar` (decision E3).
- **Tests:** the 16 test files named after the function they test (`test_lpi_ineq_sop` →
  `test_lpi_ineq`, ...).
- **Throughout:** code, help, comments, error identifiers, documentation – as in the container
  rename.
- **Kept:** names in which `sopvar` is the class (`apply_sopvar`, `compose_sopvar`, `int2sopvar`,
  `mat2sopvar`, `sopvar2sdopvar`, `settings2possopvar`, `sv_test_sopvar_*`); the PIESOS classes
  `tensopvar`, `tensopmat`; `exec_sop_settings` (it names the settings field `st.sop`); the fields
  `st.sop.*`; the test helpers `lb_full_sop`, `lb_sop_plan` (dated baselines, beside `lb_full`) and
  `bl_sop`; the local `rand_sop` (a random `sopvar`); records in `cleaning/` and `orphaned_file/`,
  dated READMEs and logs.

## 2. Inventory at 91ea628e

- `rename_table.tsv`: 71 rows – 42 toolbox functions renamed, 1 kept; 16 tests renamed, 9 kept;
  3 folders.
- 22 new names are 2025 names: 7 executives, `getController`, `getObserver`, `get_lpivar_degs`,
  `getsol_lpivar`, `lpi_eq`, `lpi_ineq`, `lpigetsol`, `lpiscript`, `lpisettings` and the six tiers.
  None of the 73 proposed names (the table and the 11 names of H2) resolves on the saved path or in
  PIETOOLS today (`namecheck_sop.m`, MATLAB R2025b, with YALMIP, MPT and SeDuMi on the path).
- Lines that name a renamed token: toolbox 851 (282 files outside `cleaning/` and `orphaned_file/`,
  266 of them .m), toolbox documents 18, tests 1,079. Kept as records: `cleaning/` 3,041,
  `orphaned_file/` 188 (`counts.txt`).
- 74 of the 266 .m files existed in 2025 (key: `@class/method`, else basename, against
  origin/PIETOOLS-2025): examples 56, demos 13, `PIETOOLS_PDE.m`, `PIETOOLS_DDE.m`,
  `PIETOOLS_DDF.m`, tests 2 (`touched_2025.txt`).
- Declan files: `PIESOS/piesos_getsol.m`.
- 103 error identifiers carry a renamed function name (`errids.txt`).
- Moves: 42 + 16 files; 3 folders and 2 test folders.

## 3. Signatures against 2025

The 20 functions with a 2025 namesake keep the 2025 order of inputs and outputs; the current
versions append optional inputs (`form`, `gain`) and a last output `info`. `lpiscript` maps each
name to the executive and form of 2025 (`'stability'` → `PIE2PDEstability` = `'Q'`,
`'l2gain-dual'` → `Hinf_gain_dual` = `'Qd'`, `'h2norm'` → `H2_norm_o` = `'o'`, ...). Differences:

| Function | 2025 | Current | A 2025 call after the rename |
|---|---|---|---|
| `getController`, `getObserver` | `(P,Z,tol)`: `inv(P,tol)`, then coefficients below `tol` dropped | `(P,Z,opts)`, method `'qp'` by default, no truncation | `tol` ignored without an error (`isfield` of a double is false) |
| `lpi_ineq` | `opts.psatz` 0/1 (add the psatz term), `opts.pure` (no R0), `opts.sep` | `opts.psatz` a term set or term codes (0 plain, 1 box product, ...); `pure`, `sep` not read (`popts` built anew, `lpi_ineq_sop.m:388`) | `pure`, `sep` ignored; a numeric `psatz` is read as codes – a valid current call |
| `lpisettings`, tiers | fields `dd1`, `options1`, ..., `sosineq_on` | `lpi.*`, `sop.*` | the executives stop with an error that names `sop.stock` (`exec_sop_settings.m:133`) |

## 4. Hazards

**H1 – 2025 inputs with another meaning.** `getController(P,Z,tol)`, `getObserver(P,Z,tol)`,
`lpi_ineq(prog,P,opts)` with `opts.pure` or `opts.sep`: accepted and ignored (decision E2). A
numeric `opts.psatz` cannot be told from current term codes by type; what 0 and 1 select against
2025 is measured in Phase 1b and documented, not caught.

**H2 – 2025 executives that became forms.** Eleven 1-D names, with the form read from the 2025 code
(the LPI in the header, the `lpiscript` map) or identified numerically (10/2026, 11/27/2026):
`PDEstability` `'P'`, `PDEstability_dual` `'Pd'`, `PIE2PDEstability` `'Q'`,
`PIE2PDEstability_dual` `'Qd'`, `Hinf_gain_coercive` `'P'`, `Hinf_gain_dual` `'Qd'`,
`Hinf_gain_dual_coercive` `'Pd'`, `H2_norm_o` `'o'`, `H2_norm_c` `'c'`, `H2_norm_o_coercive`
`'oco'`, `H2_norm_c_coercive` `'cco'` (each `PIETOOLS_` prefixed). Ten 2-D names
(`PIETOOLS_stability_2D`, `PIETOOLS_Hinf_gain_2D_non_coercive`, ...; the 2-D convention marks the
non-coercive form instead) have no checked mapping. Today each is "Undefined function" (decision
E1).

**H3 – a PIETOOLS 2025 installation on the path.** The 22 names of Section 2 join the names both
releases already share; path order decides, without a check.

**H4 – mechanics.** The MATLAB desktop PID 38352 (since 11/21) has
`executives_sopvar/PIETOOLS_H2_norm_sop.m` open: a folder rename fails with "Permission denied"
while MATLAB watches the folder (11/08), and saving that buffer would recreate the old file. Close
every MATLAB session and pause OneDrive sync before Phase 1; moves with `mv`, file by file.

**H5 – documentation.** After the pass `executives/`, `lpi_programming/`, `settings/` are the live
folders, and the retired 2025 code keeps the same names under `cleaning/private/`. Text that calls
`executives/` the retired stock path (agent file, `CLAUDE.md`, memories, READMEs) must name
`cleaning/private/executives/`.

**H6 – error identifiers.** The 103 identifiers change with their function; tests that match them
change in the same commit.

**H7 – concurrent edits.** One Declan file (`piesos_getsol.m`); other sessions share the checkout.

**H8 – records.** `cleaning/` and `orphaned_file/` are kept as written; dated baseline code is
updated so it still runs, its READMEs and logs are kept.

## 5. Decisions

Carried over from the container rename (11/26), unless the maintainer says otherwise: D1, one log
line in each of the 74 files that existed in 2025, none in files new in the release; D2, each file
that names a 2025 name of this table reviewed alone (Phase 0); D7, rename now, announce after, the
Declan file in the rename commit; the Phase 1 tools reused with this table.

Asked and answered 11/28/2026:

| | Question | Decision |
|---|---|---|
| E1 | the 2025 executives of H2 | "Error, but with specific guidance on the correct calling format for the case which led to the specific error": one file per 2025 name that stops with the current call for that name |
| E2 | the 2025 inputs of H1 | stop with an error that gives the current call (the D6 rule) |
| E3 | the folders | renamed to the 2025 names |
| E4 | commits | commit each phase when it passes its checks; push after Phase 3; announcement drafted for the maintainer |

## 6. Phases

0. Review of code that names a 2025 name of the table (callers that would bind to a renamed
   function), case by case.
1. Dry run on an export (`ren_apply.pl` with this table; checks 0, 3a–3e of the container rename),
   then the checkout: moves, folders, map.
1b. E2 (errors for the 2025 inputs of H1) and E1; the `psatz` measurement of H1.
2. Documentation: agent file, `CLAUDE.md`, READMEs, the two `sopvar` notes and their PDFs, memories.
3. Identity on exports before and after: the harness of 11/27 (55 rows, every SDP bit for bit in A,
   b, c, K, x, y; F1, F2; DEMO workspaces; class operations), batches under 20 minutes.

## 7. Effort and run time

Phase 1 edits about 1,950 lines (851 toolbox, 18 documents, 1,079 tests) plus 74 log lines; the
container rename edited 2,213. Phase 3: about 40 minutes of MATLAB per tree in batches of at most
17 minutes, plus the 3-D heat row (18 minutes per tree).
