# Announcement draft – for the maintainer to send (D7)

Drafted 11/27/2026; the maintainer sends it after the push. Recipients: Declan Jagt and the other
developers of branch `ndopvar`.

---

**Subject:** ndopvar: `copvar`, `cdopvar` renamed `opvar`, `dopvar` – commit or push local work
before pulling

Branch `ndopvar` now renames the container classes `copvar` and `cdopvar` to `opvar` and `dopvar`,
with the functions named after them (commits ae2b95c1, abfb0e83 and 2edd98e3 of 11/27/2026). `sopvar` and `sdopvar` keep
their names. The 2025 classes `opvar`, `dopvar`, `opvar2d` and `dopvar2d` have been off the path
since 11/08 (`cleaning/private/opvar/`).

**Before pulling, commit or push local work.** A merge carries committed changes into the renamed
files; an uncommitted edit to a renamed file has to be carried over by hand. Among the renamed or
edited files: 16 parser files (`mtimes_copvar` → `mtimes_opvar` in `@pde_struct/private`, and the
`InferiorClasses` of `pde_struct.m`) and the two PIESOS converters (`copvar2tensopmat`,
`tensopmat2copvar`).

Renamed (every row: `cleaning/private/rename_copvar_opvar_1125/rename_table.tsv`):

| Before | After |
|---|---|
| `copvar`, `cdopvar` (`sopvar/@copvar`, `sopvar/@cdopvar`) | `opvar`, `dopvar` (`sopvar/@opvar`, `sopvar/@dopvar`) |
| `poscopvar`, `lpivar_cdopvar`, `lpi_eq_cdopvar` | `poslpivar`, `lpivar`, `lpi_eq_dopvar` |
| `mat2copvar_sop`, `eye_copvar_sop`, `zeros_copvar_sop`, `pi4_copvar_sop` | `mat2opvar`, `eye_opvar`, `zeros_opvar`, `pi4_opvar` |
| `int2copvar`, `qp_copvar` | `int2opvar`, `qp_opvar` |
| `poscopvar_direct`, `poscopvar_lift`, `posmult_cdopvar`, `copquadvar`, `lift_copvar` | `poslpivar_direct`, `poslpivar_lift`, `posmult_dopvar`, `opquadvar`, `lift_opvar` |
| `copvar_space_list`, `parse_copvar_spaces`, `mat2copvar_grid`, `dpvar_op_copvar` and six more in `sopvar/misc` | `copvar` → `opvar` in each name |
| `copvar2tensopmat`, `tensopmat2copvar`, `mtimes_copvar` | `opvar2tensopmat`, `tensopmat2opvar`, `mtimes_opvar` |
| `poslpivar_sop` | `poslpivar_legacy` |
| `Test_copvar`, `Test_cdopvar`, `Test_qp_copvar_nd`, `rand_copvar`, `rand_cdopvar` | `Test_opvar`, `Test_dopvar`, `Test_qp_opvar_nd`, `rand_opvar`, `rand_dopvar` |

Retired, not renamed: `op2copvar_sop`, `pie2copvar` (moved to `orphaned_file/private`).

Also changed:
- Error identifiers follow the class: `copvar:...` is now `opvar:...`.
- A call in the 2025 form `poslpivar(prog,n,d,options)`, `lpivar(prog,n,d,opts)` or
  `mat2opvar(Pmat,Pdim,vars,dom)` stops with an error that gives the current call (identifiers
  `poslpivar:call2025`, `lpivar:call2025`, `mat2opvar:call2025`).
- A `.mat` file holding a `copvar` or `cdopvar` object loads as an empty double; no tracked file
  holds one.

To update local code: `cleaning/private/rename_copvar_opvar_1125/phase1/dryrun_1126/map.sed`
applies the map as whole words (GNU sed, `sed -i -f map.sed file.m`); it also rewrites comments
and header-log lines, which the commits left as written.

Verification, on exports before (4ac5c66d) and after (abfb0e83) the rename: 55 test, demo and
fixture rows give the same status; 110 SDPs – converter, LPI and executive tests, the FI library,
PIESOS, DEMO1–DEMO15 and a 3-D heat LPI with 2.7 million decision variables – are equal bit for
bit in A, b, c and K, 109 of them also in x and y; PIESIM on 73 library examples and the
operators of 172 library PIEs differ by 0
(`cleaning/private/rename_copvar_opvar_1125/phase3/README.md`).

---
