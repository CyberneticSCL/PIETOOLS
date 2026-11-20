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

`PIETOOLS_PDE_Example_Heat_Eq_with_ODE_Das.m` stays in `Stability_Examples/`, because it does not
duplicate case 8. It lacks the derivative on the PDE term of the ODE, so its ODE input is x(a), which
the BC sets to 0. Its ODE is therefore decoupled from the heat equation.
