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
