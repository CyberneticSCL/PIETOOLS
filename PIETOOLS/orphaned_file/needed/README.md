# orphaned_file/needed/

Files of the folders retired to `cleaning/private/` on 11/08/2026 (`lpi_programming/`, `opvar/`,
`executives/`, `settings/`) that code on the MATLAB path still calls (maintainer, 11/08/2026:
"anything that is needed can move into orphaned_files/needed"). This folder is on the path:
`pietools_path_update` adds it with `genpath`. Each file is a rename of a tracked file, unchanged,
under its original subpath, so `git mv orphaned_file/needed/<path> <path>` restores it. Each name
has one definition on the path.

Since 11/08/2026 the folder holds no code: each file below was here at first and left later that day.

- `settings/lpisettings.m` and the 12 presets `settings/settings_PIETOOLS_{light,heavy,veryheavy,
  stripped,custom,extreme}.m`, `settings/2D/settings_PIETOOLS_*_2D.m`: called by the container
  executives (default settings), `PIETOOLS_auto_execute_sop`, `lpiscript_sop`, `hinf_gain_1d_sop`,
  DEMO7, DEMO12, the Ch2 snippet, the cuADMM baseline suite, container tests, and the tiers of
  `settings_sopvar`. Every caller on the path now calls `lpisettings_sop`, whose tiers read copies of
  the presets they use in `settings_sopvar/private/` (with an `extreme` tier added); the only
  remaining call is in `PIETOOLS_demos/lowrank_2d_stability/tests_1d/`, which already calls retired
  functions. Now in `cleaning/private/settings/` (maintainer: "Move it now").

- `executives/utility_functions/closedLoopPIE.m`, called by `PIETOOLS_PDE.m`, DEMO12 and the Ch2,
  Ch5, Ch6 and Ch7 snippets: now in `executives_sopvar/utility_functions/`, code unchanged.
- `executives/utility_functions/get_eq_opts_2D.m`: kept for `cx_stability_eq_opts_2D` (E3), which
  only cites it in a comment and an error message; its one caller is
  `PIETOOLS_demos/lowrank_2d_stability`, which already calls retired functions. Now in
  `cleaning/private/executives/utility_functions/`; `lpi_programming_sopvar/misc/eq_opts_sopvar.m` is
  its container counterpart.

- `lpi_programming/lpi_ineq.m` and `lpigetsol.m`, for E1 (`e2e_hinf_gain_1d`),
  `test_dpvar_hinf_chain_sop` and `test_getsol_sop` (iv), which called them on a scalar `dpvar`. Those
  tests now call `lpi_ineq_sop` and `lpigetsol_sop`, which reduce to `sosineq` and `sosgetsol` exactly
  as the stock functions did, and `lpi_ineq_sop` no longer hands other classes to `lpi_ineq`; both
  files are in `cleaning/private/lpi_programming/`.
