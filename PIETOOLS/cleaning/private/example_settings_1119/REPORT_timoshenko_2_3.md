# Library cases 15 and 16: the coded Timoshenko systems are unstable

Files: `PIETOOLS_examples/Examples_Library/Stability_Examples/PIETOOLS_PDE_Ex_Timoshenko_Beam_2.m`
(case 15) and `..._Timoshenko_Beam_3.m` (case 16). Found 11/19/2026 while searching settings for every
library example: no setting certifies either row, and both systems are unstable as coded, so no
stability certificate exists. The files are not edited (git log: Sachin 2022, Declan 2024, Dan
Braghini 05/2026); this report is for their authors.

## Case 15, Timoshenko_Beam_2

- The header documents `r*II * phi_tt = E*II * phi_ss + k*aa*g * (w_s - phi) - c*phi_t`.
- The coded equation of `x3 = phi_t` is `x3_t = x2 - x4 - c*x3` in both formats: the batch matrix
  `PDE_b.A2 = [0;0;0;0]`, and the term-based `PDE_t.x{3}` has no term with `x = 4, D = 2`. The bending
  term `E*II*phi_ss` is absent. The header ends with "Unstable."
- At the default `c = 0`, the function `x(t) = (w1, t*phi1, phi1, t*phi1)`, with `phi1(s) = 2s - s^2`
  and `w1(s) = s^2 - s^3/3`, satisfies the coded PDE and its five BCs (checked symbolically in
  Wolfram), and `||x4(t)||^2 = 8t^2/15` grows without bound.
- Restoring the term (`PDE_b.A2(3) = 1`; a term `x = 4, D = 2, C = 1` in `PDE_t.x{3}`) gives the
  documented beam.

## Case 16, Timoshenko_Beam_3

- The header documents the single-equation form `0 = alp*w_ssss + bet*w_tt - gam*w_ttss + w_tttt`
  with the BCs `w(0) = w_s(0) = 0`, `w_ss(1) - w(1) = 0`, `w_sss(1) - w_s(1) = 0`, and states "stable".
- At the defaults `(alp, bet, gam) = (1, 1, 2)`, the characteristic determinant of
  `lambda^4 v + lambda^2 v - 2 lambda^2 v'' + v'''' = 0` under these BCs has the roots
  `lambda = 1.5207 ± 1.0639i, 2.2654 ± 4.3009i, 2.7484 ± 7.5350i` (FindRoot at 40 digits; Chebyshev
  collocation at 40 and 56 nodes agrees to 4 digits). The coded system grows exponentially.
- The BCs at `s = 1` need checking against the free-end conditions of the documented two-equation
  system, `phi_s(1) = 0` and `w_s(1) - phi(1) = 0`.

## Effect on the per-example settings

Both rows stay uncertified (`PHASE_C.md` §4 and §9). The infeasibility certificates the LPIs return
are correct for the coded systems.

## Literature (11/19, `PHASE_C.md` §10)

- Case 15 is the mixed hyperbolic-diffusive form of Peet, Shivakumar, Das, Weiland, IFAC-PapersOnLine
  52(2), 2019, §8.1.0.2, which has A2 = [0;0;1;0]; the paper uses c = 0 on purpose (no stability in
  these states) and reports stability for c > 0. The term was dropped from the batch format by
  72de7c6b (12/07/2022); the term format never had it outside a commented-out alternative.
- Case 16: the coded BCs are those of Shivakumar et al., IEEE TAC 70(5), 2025, §VII-B (a simulation,
  no stability claim), not the clamped-free BCs of the beam. The BCs of
  `command_line_input_Examples/stability_16_3.m` (w_sss(0) − w_tts(0) + w_s(0) = 0, w_ss(1) = w_tt(1),
  alp·w_sss(1) = gam·w_tts(1)) reproduce the spectrum of the two-equation beam, on the imaginary axis
  (Wolfram, `litcheck_1119/timoshenko_charfun.wl`): stable, not exponentially stable.

## Applied 11/19 18:05 (maintainer: "Apply all four")

- Case 15: the bending term restored in both input formats, in the library description and in
  `stability_15.m`. Certified in its dual test (`Qd`) at veryheavy with dD 0, 74,593 decision
  variables (`Qd` at stripped, light and heavy: infeasible); the primal test `Q` certifies at light
  (21,413) for c = 0, 0.1 and 1. The file returns this st and names its form.
- Case 16: unchanged (its BCs come from a published simulation; the choice between the beam BCs of
  `stability_16_3.m` and a corrected header is open).

## Case 16 removed, 11/19 18:40 (maintainer: "remove it from the stability examples")

The file, its GUI file and the command-line copies with the same BCs (`stability_16.m`,
`stability_16_2.m`) are in `orphaned_file/private/`; the library entry stops with an error naming the
removal. The beam BCs remain in `command_line_input_Examples/stability_16_3.m`. Primary source read
(arXiv 1812.06794v1, §8.1.0.2): the case 15 form with A2 = [0;0;1;0], not certified there at c = 0,
stable for c > 0.
