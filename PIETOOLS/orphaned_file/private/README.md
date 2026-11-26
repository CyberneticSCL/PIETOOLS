# orphaned_file/private/

Code that no toolbox function calls, kept for reference off the MATLAB path: `genpath`
(`pietools_path_update`) skips a private folder, and no function in `orphaned_file/` resolves into it.

## Moved 11/19/2026

Moved with `mv`, so the git index still lists the original paths; `mv` back restores them.

- **Library case 16**, from `PIETOOLS_examples/Examples_Library/Stability_Examples/`:
  `PIETOOLS_PDE_Ex_Timoshenko_Beam_3.m` and `PIETOOLS_PDE_Ex_Timoshenko_Beam_3_GUI.mat`. Its
  command-line copies with the same BCs, `stability_16.m` and `stability_16_2.m`, came from
  `.../command_line_input_Examples/`.
  - Maintainer: "remove it from the stability examples".
  - The example is the one-equation Timoshenko beam
    w_tttt + w_tt − 2 w_ttss + w_ssss = 0 with BCs w(0) = w_s(0) = 0, w_ss(1) = w(1),
    w_sss(1) = w_s(1). These are the BCs of a simulation in Shivakumar, Das, Weiland, Peet, IEEE TAC
    70(5), 2025, Sec. VII-B.
  - The system has the eigenvalues 1.5207 ± 1.0639i, 2.2654 ± 4.3009i and 2.7484 ± 7.5350i, so it
    is unstable, and its header claim "stable" is false.
  - The library entry now stops with an error that names the removal.
  - The beam BCs are in `stability_16_3.m`, which stays in place.
- **Predecessors that duplicate library examples**, from `Stability_Examples/` (maintainer: "move
  them to orphaned code if duplicates"):
  - `PIETOOLS_PDE_Example_Saba.m`, with `PIETOOLS_PDE_Example_Saba_Ex_1_GUI.mat` to `_Ex_5_GUI.mat`:
    library case 4. Its case 4.5 keeps the typo sig2 = 2.048, corrected to 1.048 in the library
    file.
  - `PIETOOLS_PDE_Example_Lamare.m`, with `PIETOOLS_PDE_Example_Lamare_Ex_5_1_GUI.mat` to
    `_Ex_5_4_GUI.mat`: library case 2.
  - `PIETOOLS_PDE_Example_Transport_Eq.m`: library case 1 mirrored, x_t = x_s with the BC at s = 1.
    Its GUI file stays in `Stability_Examples/`, because `stability_1.m` loads it.

## Moved 11/21/2026

- **Predecessor of library case 8**, from `Stability_Examples/` (maintainer, 11/21/2026):
  `PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m`.
  - Its batch format is case 8 with the sign +Bxr of the paper, corrected to −Bxr in the library
    file; with +Bxr the system is unstable at its default lam = pi^2 − 1.
  - Its term-based format lacks the derivative on the PDE term of the ODE, so its ODE input is x(a),
    which the BC sets to 0, and its ODE is decoupled from the heat equation.
  - With GUI = 1 it loads `PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das_GUI.mat`, which does not exist.

## Moved 11/22/2026

- **Periodic-BC map**, from `converters/internal/` (maintainer, 11/22/2026: "move to orphaned code"):
  `compute_PIETOOLS_PDE_periodic_map.m` (DB, 10/28/2026).
  - It builds its three operators as `opvar` objects (`opvar()` for Fop, Tx and Kop). That class has
    been retired since 11/08/2026, so the file stops at its first `opvar()` call.
  - Its only caller is `converters/legacy/convert_PIETOOLS_PDE_2D.m` (line 502). That file belongs to a
    legacy 2-D converter group that nothing outside the group calls (audit of 11/22/2026,
    `cleaning/private/dependency_audit_1122/AUDIT.md`, §2).
- **Converters from the opvar classes**, from `converters/operator_converters/` (maintainer,
  11/22/2026: move the branches that run only for opvar inputs to orphaned code): `opvar2copvar.m`
  (MMP, 09/19/2026), `opvar2sopvar.m` and `opvar2d2sopvar.m` (DJ, 05/27/2026), `dopvar2sdopvar.m`
  (DJ, 10/05/2026).
  - Each converts an object of a class retired on 11/08/2026.
  - Their callers were the branches below and each other. The one call left is in the legacy 2-D
    converter group (`convert_PIETOOLS_PDE_nonlinear_2D.m`, three `opvar2sopvar` calls).
- **Branches for opvar inputs**, removed from seven live files on the same instruction. The files as
  they were before the removal are in `opvar_branches_1122/` here, with the path below `PIETOOLS/`
  (`@X` written `X_class`, so that no copy can define a class):
  - `converters/operator_converters/op2copvar_sop.m`: the conversion of opvar, dopvar and opvar2d. It
    now returns a container as is and errors for any other class.
  - `parser/piess.m`: the call of `piess_opvar` when the operators are opvar objects and none is a
    copvar.
  - `PIESOS/piesos_getsol.m`: the branches for dopvar, dopvar2d (`getsol_lpivar`) and ndopvar, and the
    opvar, opvar2d and nopvar terms of the pass-through condition.
  - `PIESOS/piesos_ineq.m`: the branch that sent an opvar-class F to `lpi_ineq`.
  - `PIESOS/polyopvar/@tensopmat/tensopmat.m` (two branches), `@tensopmat/initialize.m` and
    `@tensopvar/initialize.m` (one each): the conversion of opvar and opvar2d blocks by
    `opvar2sopvar`.
  - Checked on 11/22/2026 by running DEMO1, DEMO2, DEMO7, the Ch5 snippet, the three PIESOS demos
    and six tensopmat and tensopvar tests on HEAD and on the edited tree. All 13 run, and their output
    is the same apart from solver times.
- **Legacy 2-D converter**, into `converters/legacy/` here, with the same paths (maintainer,
  11/22/2026: "move the legacy 2-D converters to orphaned code"): `convert_PIETOOLS_PDE_2D.m`,
  `convert_PIETOOLS_PDE_nonlinear_2D.m`, `pie_pde_maps_2D.m`, and from `2D/`
  `Compute_Tmap_PDE_2D_Separable.m`, `PIETOOLS_FTC_Expansion_2D.m`, `expand_PIETOOLS_PDE_BCs.m` (DJ,
  2022–2026) and `compute_PIETOOLS_PDE_periodic_map_2D.m` (DB, 10/29/2026).
  - The seven files called only each other. They build `opvar2d` operators (`opvar2d()`,
    `clean_opvar`, `mat2opvar`) and call `opvar2sopvar` and `compute_PIETOOLS_PDE_periodic_map`,
    which are also here now. 2-D PDEs convert through `convert` on `pde_struct` (DJ, since 11/02/2026).
  - Checked on 11/22/2026: three 2-D library examples (`2D_Heat_Eq_with_ODE`, `2D_KISS_Model`,
    `2D_NS_Antonelli`) convert to `copvar` PIEs before and after. The PIESOS test
    `Test_convert_PDE_nonlinear_varnames` called the nonlinear converter directly and stopped on the
    retired `opvar`. Since 146cbb70 (DJ, 11/24/2026) it calls `convert`, and no file outside
    `cleaning/` and `orphaned_file/` calls the seven files.
- **Low-rank 2-D certify family**, into `PIETOOLS_demos/lowrank_2d_stability/` here, with the same
  path (maintainer, 11/22/2026: "port or retire"; retired by the low-rank session, CC):
  `pielr_certify.m`, `pielr_certify_pos.m`, `pielr_tensor.m`, `build_stab_2d_st2.m`, `face_n1.mat`,
  `demo_lowrank_2d.m`, `demo_your_system.m`, `demo_poincare_2d.m`, and from `private/`
  `opcheck_2d.m`, `opcheck_pos2d.m`, `pielr_evalgram.m`, `restrict_solve_pos.m`, `sed_ref_pos.m`,
  `nb_2d.m`, `nb_monoms.m`.
  - They build and check `opvar2d` programs (`poslpivar_2d`, `lpi_eq_2d`, `lpigetsol`,
    `clean_opvar`, `get_eq_opts_2D`), which left the path on 11/08/2026, so each stopped at its first
    call to them (audit §2).
  - What they added over `pielr_solve`: a stored face (`face_n1.mat`, in `opvar2d` Gram coordinates)
    lifted to n states by `pielr_tensor`, and the low-rank certificate of a given 2-D positive operator
    (`pielr_certify_pos`). Face transfer was abandoned on 09/20/2026. A port would need the face
    regenerated in the container Gram layout and a container gate for the positive-operator case.
  - 2-D stability stays available through `pielr_solve` (container builders). The live files edited
    with this move: `pielr_path.m` (assert list), `private/restrict_solve.m` (the adapter is
    required), `private/pielr_evalop.m` (two unreachable `opvar` branches deleted).

## Moved 11/26/2026

Moved with `mv`, paths below `PIETOOLS/` kept, ahead of the rename of `copvar` and `cdopvar` to
`opvar` and `dopvar` (maintainer decisions C3, C4 and D4; plan and review in
`cleaning/private/rename_copvar_opvar_1125/`).

- **Legacy terms converter**, into `converters/legacy/` here: `convert_PIETOOLS_PDE_terms_internal.m`
  (SS, 12/05/2022).
  - No file calls it. It declares its operators by lines such as `opvar Aop_PDE Bop_PDE;`, so since
    the class left the path on 11/08/2026 it stopped at the first of them. After the rename those
    lines would build empty containers and fail later, at the `.dim` assignments (measured).
- **Legacy state/terms input format**, into `parser/legacy/` here: the class folders `@state` (20
  files) and `@terms` (14 files).
  - No live code constructs a `state` or `terms` object, and `terms2pde_struct` raises "Objects of
    type 'terms' (a legacy input format) are no longer supported" (DJ, 11/25/2026). Since 11/08/2026
    the first operation on a state, `x+x` or `-x`, stopped at `opvar T;` (measured).
  - After the rename those lines would have run with a warning, and a legacy script would have got
    further, on a path never tested.
  - `parser/legacy/@sys` stays, since `convert_PIETOOLS_PDE` accepts a `sys` object. The tests
    `isa(x,'state')` and `isa(x,'terms')` in `@sys` and in the `pde_struct` operators now return
    false.
- **Converters that returned containers unchanged**, into `converters/operator_converters/` here:
  `op2copvar_sop.m` (its input unchanged since 11/22/2026) and `pie2copvar.m` (a container PIE
  unchanged; otherwise `initialize`, which accepts only containers since b1211608).
  - Renamed, they would read `op2opvar` and `pie2opvar`; `pie2opvar` named the opposite conversion
    until 6757b6d3.
  - Their calls were deleted from DEMO7, DEMO12 and the three copies of the baseline harness
    `hd_qform.m`. `pie2copvar` left the path assertions of the four `lb_setup.m` copies, and
    `op2copvar_sop` joined the retired names in `p0_path.m`.
