# orphaned_file/needed/

Files of the folders retired to `cleaning/private/` on 11/08/2026 (`lpi_programming/`, `opvar/`,
`executives/`, `settings/`) that code on the MATLAB path still calls (maintainer, 11/08/2026:
"anything that is needed can move into orphaned_files/needed"). This folder is on the path:
`pietools_path_update` adds it with `genpath`. Each file is a rename of a tracked file, unchanged,
under its original subpath, so `git mv orphaned_file/needed/<path> <path>` restores it. Each name
has one definition on the path.

| File | Called by |
|---|---|
| `settings/lpisettings.m` | the container executives `PIETOOLS_*_sop` (default settings), `PIETOOLS_auto_execute_sop`, `lpiscript_sop`, `lpi_programming_sopvar/examples/hinf_gain_1d_sop.m`, DEMO7, DEMO12, the Ch2 snippet, the cuADMM baseline suite, container tests |
| `settings/settings_PIETOOLS_{light,heavy,veryheavy,stripped,custom,extreme}.m` and `settings/2D/settings_PIETOOLS_*_2D.m` (12) | `lpisettings`, one 1-D and one 2-D file per preset; `settings_sopvar/settings_sopvar_{light,heavy,veryheavy,stripped,custom}.m` |
| `executives/utility_functions/closedLoopPIE.m` | `PIETOOLS_PDE.m`, DEMO12, the Ch2, Ch5, Ch6 and Ch7 snippets |
| `executives/utility_functions/get_eq_opts_2D.m` | `cx_stability_eq_opts_2D`, so E3 (`e2e_stability_2d`) |

None of them calls a retired function: `closedLoopPIE` only tests for the opvar classes by name
(`isa`).

`lpi_programming/lpi_ineq.m` and `lpigetsol.m` were here until later on 11/08/2026, for E1
(`e2e_hinf_gain_1d`), `test_dpvar_hinf_chain_sop` and `test_getsol_sop` (iv), which called them on a
scalar `dpvar`. Those tests now call `lpi_ineq_sop` and `lpigetsol_sop`, which reduce to `sosineq` and
`sosgetsol` exactly as the stock functions did, and `lpi_ineq_sop` no longer hands other classes to
`lpi_ineq`; both files are in `cleaning/private/lpi_programming/`.
