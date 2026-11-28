# Announcement draft – for the maintainer to send (E4)

Drafted 11/28/2026; the maintainer sends it after the push. Recipients: Declan Jagt and the other
developers of branch `ndopvar`.

---

**Subject:** ndopvar: the `_sop` suffix dropped, `executives_sopvar` is `executives` – commit or
push local work before pulling

Branch `ndopvar` now drops the suffix `_sop` from the executives and from the LPI and settings
functions, taking the PIETOOLS 2025 name where one exists, and renames the folders
`executives_sopvar`, `lpi_programming_sopvar` and `settings_sopvar` to `executives`,
`lpi_programming` and `settings` (commits fe6d3530, 0d02ac08 and 182a2c49 of 11/28/2026). The
retired 2025 code keeps those folder names under `cleaning/private/`, off the path.

**Before pulling, commit or push local work.** A merge carries committed changes into the moved
files; an uncommitted edit to a moved file has to be carried over by hand. Among the edited files:
`PIESOS/piesos_getsol.m` (a call of `getsol_lpivar`).

Renamed (every row: `cleaning/private/rename_sop_1128/rename_table.tsv`):

| Before | After |
|---|---|
| `PIETOOLS_stability_sop`, `PIETOOLS_Hinf_gain_sop`, `PIETOOLS_H2_norm_sop`, ... (10 executives) | `PIETOOLS_stability`, `PIETOOLS_Hinf_gain`, `PIETOOLS_H2_norm`, ... |
| `lpiscript_sop`, `lpi_ineq_sop`, `lpi_eq_sop`, `lpigetsol_sop`, `getsol_lpivar_sop` | `lpiscript`, `lpi_ineq`, `lpi_eq`, `lpigetsol`, `getsol_lpivar` |
| `lpisettings_sop`, `lpisettings_stock_sop`, `settings_sopvar_<tier>` | `lpisettings`, `lpisettings_stock`, `settings_PIETOOLS_<tier>` |
| `getController_sop`, `getObserver_sop`, `get_lpivar_degs_sop`, `pie_witness_sop`, `pie_disc_sop` | `getController`, `getObserver`, `get_lpivar_degs`, `pie_witness`, `pie_disc` |
| `exec_info_sop`, `exec_pie_sop`, `exec_tools_sop`, `lpi_solve_loop_sop`, `poslpivar_settings_sop`, ... | the same without `_sop` |
| `eq_opts_sopvar` | `get_eq_opts` |
| `test_executives_sop`, `test_lpi_ineq_sop`, `test_call2025_sop`, ... (16 tests) | the same without `_sop` |

Kept: `exec_sop_settings` (it names the settings field `st.sop`), the names of the `sopvar` class
(`apply_sopvar`, `int2sopvar`, ...), and the helpers of the dated baselines.

Also changed:
- Error identifiers follow their function (`lpiscript_sop:...` is now `lpiscript:...`, 103 in all).
- The 21 executives of PIETOOLS 2025 whose LPI is now a form (`PIETOOLS_PIE2PDEstability`,
  `PIETOOLS_Hinf_gain_dual`, `PIETOOLS_H2_norm_c_coercive`, `PIETOOLS_stability_2D`, ...) stop with
  an error that gives the call of this release, e.g. `PIETOOLS_stability(PIE,settings,'Q')`
  (`executives/retired_2025/`).
- `getController(P,Z,tol)` and `getObserver(P,Z,tol)` with a numeric third input, and `lpi_ineq`
  with `opts.pure` or `opts.sep`, stop with an error that gives the current call. A numeric
  `opts.psatz` of `lpi_ineq` is read as term codes: the 2025 `psatz = 1` is `'auto'` here.

To update local code: `cleaning/private/rename_sop_1128/phase1/dryrun_1128/map.sed` applies the map
as whole words (GNU sed, `sed -i -f map.sed file.m`); it also rewrites comments and header-log
lines, which the commits left as written.

Verification, on an export of the renamed tree against the run of 11/27 on the same code before
the rename: 56 test, demo and fixture rows give the same status; 110 SDPs – converter, LPI and
executive tests, the FI library, PIESOS, DEMO1–DEMO15 and a 3-D heat LPI with 2.7 million
decision variables – are equal bit for bit in A, b, c and K, and in x and y wherever captured;
PIESIM on 73 library examples and the operators of 172 library PIEs differ by 0
(`cleaning/private/rename_sop_1128/phase3/README.md`).

---
