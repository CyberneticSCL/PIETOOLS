# G4_examples_a – ledger of the comment pass of 11/21/2026

Agent G4_examples_a. List: `G4_examples_a.tsv` (37 files, all class `condense`). Procedure of `BRIEF.md`:
prep (private copy writing its diff files to the scratch folder `G4/`, logic unchanged), header note by hand,
`cpass.pl apply`, hand edits listed per file, `cpass.pl check`, `codeident.sh`, `git ls-files --eol`
(all files `w/crlf` before and after), `git diff --numstat` = `git diff -w --numstat`.

Apply counts: dropped (was) lines / kept (was) lines / retagged markers.

## Summary

- Classes: 37 condense, edited; 0 asnew, 0 review, 0 deferred; no file skipped (all clean in
  `git status` before the first edit). TAG `MMP, 11/21/2026` for 35 files; `DJ, MMP, 11/21/2026` for
  `2D/PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay.m`; `DB, MMP, DJ, 11/21/2026` for
  `examples_PDE_library_PIETOOLS.m`.
- check PASS on 34 files. Deliberate failures, all tool limitations: `examples_DDE_library_PIETOOLS.m` (4)
  and `examples_NDSDDF_library_PIETOOLS.m` (3), (was) lines kept because they are 2025 commented example
  code; `examples_PDE_library_PIETOOLS.m` (5), the 2025 log entries reworded on 11/03/2026 whose 2025 date
  formats the tool does not match. codeident OK 37/37; eol `w/crlf` 37/37 as before; numstat = numstat -w
  37/37. MATLAB parse check (`syntax_check`, lock runner, 06:19): 37 files, 0 with parse problems.
- Rerun with the updated `cpass.pl` (06:14 version: `DJ (CC)` initials, markers without a blank, per-process
  diff files), on request of the coordinator: no file of the list has `(CC)` or a marker without a blank
  (current files and HEAD); fresh prep gives the same code diff as the prep used for every note (37/37;
  the earlier prep ran from a private copy with private diff files, so no cross-agent mix-up); apply
  rerun on 35 files leaves them byte-identical; on the DDE and NDSDDF libraries apply was previewed on a
  copy only (`CPASS_WF`): it would drop exactly the 4 and 3 hand-restored (was) lines. No note changed.
- Open questions for the maintainer: (1) keep the 7 (was) lines of 2025 commented example code in the DDE
  and NDSDDF libraries (kept) or drop them as `apply` would; (2) `PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay.m`:
  the sign change of DJ carries no marker or (was) line (none added); (3) `examples_PDE_library_PIETOOLS.m`:
  the case-16 error message (code) contains the date 11/19/2026; line 1 marker removed and nested
  example-72 BEGIN/END pair removed (development history) – confirm.

## Files

1. `PIETOOLS_examples/Examples_Library/2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDD.m` –
   condense, TAG `MMP, 11/21/2026`. Note:
   `% MMP, 11/21/2026: Changes since PIETOOLS 2025: second output st, the` /
   `%   recommended settings of the H-infinity gain test; Hinf_gain_form = 'P'` /
   `%   set in the base workspace; help revised, license added.`
   Apply 0 / 1 / 9 (= prep). Internal notes removed: the two header entries (search of 11/19, maintainer)
   and the release-pass entry. check PASS, codeident OK, eol w/crlf, numstat OK.
2. `.../2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DDDN.m` – condense, TAG
   `MMP, 11/21/2026`. Note: same wording as 1. Apply 0 / 1 / 8 (= prep). Removed: the same two entries.
   PASS, OK, w/crlf, numstat OK.
3. `.../2D/Hinf_Gain_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN.m` – condense, TAG
   `MMP, 11/21/2026`. Note: as 1, with `help revised, its L2-gain formula corrected` /
   `%   to mu_{m,n}; license added.` (the help formula now matches the code). Apply 0 / 1 / 9 (= prep).
   PASS, OK, w/crlf, numstat OK.
4. `.../2D/Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_2D_ReactionDiffusion_DNDN_EdgeObserve.m` – condense,
   TAG `MMP, 11/21/2026`. Note: `second output st, the recommended settings of the H-infinity optimal`
   `estimator synthesis; help revised, its outputs y1, y2 corrected to the integrals the code takes;`
   `license added.` Apply 0 / 1 / 5 (= prep). PASS, OK, w/crlf, numstat OK.
5. `.../2D/PIETOOLS_PDDE_Ex.m` – condense, TAG `MMP, 11/21/2026`. Note: `second output st, the recommended`
   `settings of the H2 norm test; H2_norm_form = 'oco' set in the base workspace; the delays tau are set`
   `before params is evaluated, so params can change them; the unused pvar s t removed; help corrected`
   `(states x2, x3, output over [0,pi]), license added.` Apply 2 / 3 / 10 (= prep; the 2 dropped (was)
   lines are moved comments of 11/20). Removed: entries of 11/19 (search, maintainer), 11/20, release pass.
   PASS, OK, w/crlf, numstat OK.
6. `.../2D/PIETOOLS_PDE_Ex_2D_Heat_Eq_with_ODE.m` – condense, TAG `MMP, 11/21/2026`. Note: `optional second`
   `output st, the recommended settings of the stability test; stability_form = 'Q' set in the base`
   `workspace; help revised (OUTPUT lists no PDE_b, which is not returned), license added.`
   Apply 0 / 1 / 9 (= prep). Removed: five entries of 11/18-11/19 (tests, maintainer, settings notes).
   PASS, OK, w/crlf, numstat OK.
7. `.../2D/PIETOOLS_PDE_Ex_2D_KISS_Model.m` – condense, TAG `MMP, 11/21/2026`. Note: `second output st, the`
   `recommended settings of the stability test; stability_form = 'Q' set in the base workspace; help`
   `revised (OUTPUT lists no PDE_b, which is not returned), license added.` Apply 0 / 1 / 7 (= prep).
   PASS, OK, w/crlf, numstat OK.
8. `.../2D/PIETOOLS_PDE_Ex_2D_NS_Antonelli.m` – condense, TAG `MMP, 11/21/2026`. Note: same wording as 7.
   Apply 0 / 1 / 6 (= prep). Removed: entry of 11/19 ("no measured setting certifies it", maintainer).
   PASS, OK, w/crlf, numstat OK.
9. `.../2D/PIETOOLS_PDE_Ex_2D_Parabolic_Eq.m` – condense, TAG `MMP, 11/21/2026`. Note: as 6 (optional second
   output st, stability_form = 'Q', help revised (no PDE_b), license added). Apply 0 / 1 / 10 (= prep).
   Removed: five entries of 11/18-11/19. PASS, OK, w/crlf, numstat OK.
10. `.../2D/PIETOOLS_PDE_Ex_2D_Reaction_Diffusion_Eq.m` – condense, TAG `MMP, 11/21/2026`. Note: as 9 with
   `stability_form = 'P'` and `help revised (OUTPUT lists no PDE_b, which is not returned; stability`
   `statement added), license added.` Apply 2 / 1 / 10 (= prep; dropped: st lines of 11/18 kept as (was)
   on 11/19). PASS, OK, w/crlf, numstat OK.
11. `.../2D/PIETOOLS_PDE_Ex_2D_Telegraph_Eq.m` – condense, TAG `MMP, 11/21/2026`. Note: as 7, with `help revised`
   `(OUTPUT lists no PDE_b, which is not returned; coefficient C^2/(2*lam) corrected), license added.`
   Apply 0 / 1 / 6 (= prep). Removed: entry of 11/19 ("no measured setting certifies it", maintainer).
   PASS, OK, w/crlf, numstat OK.
12. `.../2D/PIETOOLS_PDE_Ex_2D_Transport_Eq.m` – condense, TAG `MMP, 11/21/2026`. Note: as 9 (optional st,
   stability_form = 'Q', help revised (no PDE_b), license added). Apply 0 / 1 / 10 (= prep). PASS, OK,
   w/crlf, numstat OK.
13. `.../2D/PIETOOLS_PDE_Ex_2D_Wave_Eq.m` – condense, TAG `MMP, 11/21/2026`. Note: as 12, with `(OUTPUT lists`
   `no PDE_b, which is not returned; stability statement added)`. Apply 3 / 1 / 9 (= prep; dropped: st lines
   of 11/18 kept as (was) on 11/19). PASS, OK, w/crlf, numstat OK.
14. `.../2D/PIETOOLS_PDE_Ex_Coupled_ODE_Heat_Eq_w_Delay.m` – condense, TAG `MMP, 11/21/2026`. Note: same wording
   as 7. Apply 0 / 1 / 8 (= prep). PASS, OK, w/crlf, numstat OK.
15. `.../2D/PIETOOLS_PDE_Ex_Heat_Eq_w_Delayed_Boundary_Input.m` – condense, TAG `MMP, 11/21/2026`. Note: as 7,
   with `(OUTPUT lists no PDE_b, which is not returned; stability statement added)`. Apply 0 / 1 / 6
   (= prep). Removed: entry of 11/19 ("no measured setting certifies it", maintainer). PASS, OK, w/crlf,
   numstat OK.
16. `.../2D/PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay.m` – condense, TAG `DJ, MMP, 11/21/2026` (DJ 11/12/2026
   commit 2ccb8943, sign change without a note; MMP from 11/17). Note:
   `% DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: sign of the delayed` /
   `%   reaction term corrected to -a1*x(t-tau,s1), the PDE of the help;` /
   `%   optional second output st, the recommended settings of the stability` /
   `%   test; stability_form = 'Q' set in the base workspace; help revised` /
   `%   (OUTPUT lists no PDE_b, which is not returned), license added.`
   Apply 4 / 1 / 12 (= prep; dropped: comment and st lines of 11/17-11/18 kept as (was)). Removed: seven
   entries of 11/17-11/19 (tests, maintainer, settings history). PASS, OK, w/crlf, numstat OK.
   Open: the DJ code line (`PDE_t = [... -a1*subs(x2,s2,tau);`) has no marker and no (was) line; none added.
17. `.../2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay.m` – condense, TAG `MMP, 11/21/2026`. Note: same wording
   as 7. Apply 0 / 1 / 6 (= prep). Removed: entry of 11/19 ("no measured setting certifies it",
   maintainer). PASS, OK, w/crlf, numstat OK.
18. `.../2D/PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay_2D.m` – condense, TAG `MMP, 11/21/2026`. Note: as 7, with
   `(OUTPUT lists no PDE_b, which is not returned; the BCs assume b = 1)`. Apply 0 / 1 / 6 (= prep).
   PASS, OK, w/crlf, numstat OK.
19. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance.m` – condense, TAG
   `MMP, 11/21/2026`. Note: `third output st, the recommended settings of the H-infinity gain test;`
   `Hinf_gain_form = 'Q' set in the base workspace; help revised, license added.` Apply 0 / 1 / 7 (= prep).
   PASS, OK, w/crlf, numstat OK.
20. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Diffusive_Eq_with_Distributed_Disturbance_2.m` – condense, TAG
   `MMP, 11/21/2026`. Note: same wording as 19. Apply 0 / 1 / 7 (= prep). PASS, OK, w/crlf, numstat OK.
21. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_Distributed_Disturbance.m` – condense, TAG
   `MMP, 11/21/2026`. Note: `third output st, the recommended settings of the H-infinity gain test;`
   `Hinf_gain_form = 'P' set in the base workspace; help revised: the stated gain 0.3333, that of`
   `x_t = x_ss + w, corrected to gamma = 5/24 = 0.2083, the gain at frequency 0 with the disturbance s*w`
   `(steady state x = s/2 - s^3/6); license added.` Apply 1 / 1 / 10 (= prep; dropped: the (was) help line
   `gamma = 0.3333`, a comment, its value now in the note). Removed: the LPI value 0.20847 (measurement).
   PASS, OK, w/crlf, numstat OK.
22. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Parabolic_Eq_with_Disturbance.m` – condense, TAG `MMP, 11/21/2026`.
   Note: `third output st, ...; Hinf_gain_form = 'P' set in the base workspace; the warning that SeDuMi`
   `requires the veryheavy settings removed; help revised (stated gain 15.147 corrected to 23.706,`
   `attained at frequency 0), license added.` Apply 0 / 2 / 7 (= prep; the 2025 disp line kept as (was)).
   Removed: the reason "SeDuMi certifies the gain at the settings st" (a test result). PASS, OK, w/crlf,
   numstat OK.
23. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Disturbance.m` – condense, TAG
   `MMP, 11/21/2026`. Note: `third output st, ...; Hinf_gain_form = 'Q' set in the base workspace; the`
   `SeDuMi warning states that the settings st do not certify the gain with SeDuMi and the light tier`
   `does; help revised (stated gain 8.1069 corrected to 8.2141, with its formula), license added.`
   Apply 0 / 2 / 9 (= prep). Removed: "tuned for MOSEK", the tiers tested. PASS, OK, w/crlf, numstat OK.
24. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance.m` – condense, TAG `MMP, 11/21/2026`.
   Note: as 19 with `Hinf_gain_form = 'P'`. Apply 0 / 1 / 8 (= prep). PASS, OK, w/crlf, numstat OK.
25. `.../Hinf_Gain_Examples/PIETOOLS_PDE_Ex_Wave_Eq_Tip_Damped.m` – condense, TAG `MMP, 11/21/2026`. Note: as
   24, with `help revised (stated gain for k = 0.5 corrected from 2 to 0.8419)`. Apply 0 / 1 / 7 (= prep).
   PASS, OK, w/crlf, numstat OK.
26. `.../Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_with_Control.m` – condense, TAG
   `MMP, 11/21/2026`. Note: `second output st, the recommended settings of the H-infinity optimal`
   `controller synthesis; help revised (regulated output z corrected to the one the code sets), license`
   `added.` Apply 0 / 1 / 11 (= prep). PASS, OK, w/crlf, numstat OK.
27. `.../Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Boundary_ODE_Control.m` – condense,
   TAG `MMP, 11/21/2026`. Note: as 26 with `help revised, license added.` Apply 0 / 1 / 9 (= prep). PASS, OK,
   w/crlf, numstat OK.
28. `.../Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Reaction_Diffusion_Eq_with_Controlled_Input.m` –
   condense, TAG `MMP, 11/21/2026`. Note: `third output st, ...; help revised (title of the system`
   `corrected), license added.` Apply 0 / 1 / 9 (= prep). PASS, OK, w/crlf, numstat OK.
29. `.../Hinf_Optimal_Controller_Examples/PIETOOLS_PDE_Ex_Wave_Eq_Boundary_ODE_Control.m` – condense, TAG
   `MMP, 11/21/2026`. Note: as 26 with `help revised (order of the regulated outputs corrected to the`
   `code), license added.` Apply 0 / 1 / 8 (= prep). PASS, OK, w/crlf, numstat OK.
30. `.../Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_Heat_Eq_with_Observed_Output.m` – condense, TAG
   `MMP, 11/21/2026`. Note: `third output st, the recommended settings of the H-infinity optimal estimator`
   `synthesis; help revised, license added.` Apply 0 / 1 / 9 (= prep). PASS, OK, w/crlf, numstat OK.
31. `.../Hinf_Optimal_Observer_Examples/PIETOOLS_PDE_Ex_Transport_Eq_with_Observed_Output.m` – condense, TAG
   `MMP, 11/21/2026`. Note: as 30 with `help revised (disturbance term corrected to the w(t) of the code)`.
   Apply 0 / 1 / 8 (= prep). PASS, OK, w/crlf, numstat OK.
32. `.../Examples_Library/PIETOOLS_PDE_Ex_Transport_Eq_with_Disturbance2.m` – condense, TAG `MMP, 11/21/2026`.
   Note: `second output st, the recommended settings of the H2 norm test; H2_norm_form = 'o' set in the`
   `base workspace; the unused pvar s1 s2 removed (the PDE is in s); help revised, license added.`
   Apply 0 / 2 / 6 (= prep; `%pvar s1 s2` kept as (was)). PASS, OK, w/crlf, numstat OK.
33. `.../Examples_Library/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion.m` – condense (66 % of 2025 code left), TAG
   `MMP, 11/21/2026`. Note: `second output st, the recommended settings of the H2 optimal estimator`
   `synthesis; a nonzero GUI raises an error, since no GUI file exists for this example (the 2025 code`
   `loaded the GUI file of the Euler-Bernoulli beam example); help revised (OUTPUT lists no PDE_b, which`
   `is not returned; no parameter can be set), license added.` Apply 3 / 7 / 11 (= prep; kept: the 2025 GUI
   block, 6 lines, and the signature; dropped: deg_Q line and comment text of 11/20). Removed: entries of
   11/20 (deg_Q history) and 11/21 (Initial coding correction note). PASS, OK, w/crlf, numstat OK.
34. `.../Examples_Library/PIETOOLS_PDE_Ex_Unstable_ReactionDiffusion_Control.m` – condense, TAG
   `MMP, 11/21/2026`. Note: `second output st, the recommended settings of the H2 optimal controller`
   `synthesis; help revised (control term s*u added to the PDE, as in the code), license added.`
   Apply 3 / 1 / 10 (= prep). Removed: entry of 11/20 (deg_Q history). PASS, OK, w/crlf, numstat OK.
35. `PIETOOLS_examples/examples_DDE_library_PIETOOLS.m` – condense, TAG `MMP, 11/21/2026`. Note (after the 2025
   lines `Initial coding MP - 10_01_2020`, `updated to DDE data structure 5_3_2021`):
   `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the active block (Example 2b` /
   `%   of [13], H2 norm) sets st, its recommended settings; the other blocks` /
   `%   carry, commented, their st and the form of their test (stability_form,` /
   `%   Hinf_gain_form); the wind tunnel block sets z = [x1; x2; .1u] through` /
   `%   DDE.C1, C1i, D11, D12 (the fields C0, C{1}, D1, D2 it set are not read,` /
   `%   so z was 0) and its delay as the discrete delay Ai{1} of [5] (was` /
   `%   Adi{1}); help and license added.`
   Apply 18 / 0 / 130 (= prep). Hand edit after apply: the 4 (was) lines of the wind tunnel block
   (`% %  DDE.Adi{1}=...`, `% %  DDE.C0=...; DDE.C{1}=...`, `% %  DDE.D1=...`, `% %  DDE.D2=...`) restored with
   `% MMP, 11/21/2026 (was)`: they are the 2025 lines of that commented example block (2025 code of the
   example, one comment level more), so the rule "keep commented code of the prior release" applies; apply
   drops them because the 2025 lines are themselves comments. Net: 14 dropped (the st lines of 11/18 kept
   as (was) on 11/19), 4 kept. Removed: entries of 11/18-11/19 (measurements, cleaning/ paths, D100/D115,
   maintainer). check: 4 failures, exactly these 4 lines ("(was) code not in 2025") – tool limitation,
   deliberate; codeident OK, w/crlf, numstat OK.
36. `PIETOOLS_examples/examples_NDSDDF_library_PIETOOLS.m` – condense, TAG `MMP, 11/21/2026`. Note (after the 2025
   line `Initial coding MP - 10_01_2020`):
   `% MMP, 11/21/2026: Changes since PIETOOLS 2025: the active block (multiple` /
   `%   showering people, DDF, H-infinity control) sets the regulated output of` /
   `%   option 2, z = [sum_i T1i; .1*sum_i ui] (no output was set), the delays` /
   `%   tau(i) = i (were n), and st, its recommended settings; the other blocks` /
   `%   carry, commented, their st and the form of their test (stability_form);` /
   `%   help and license added.`
   Apply 3 / 1 / 75 (= prep). Hand edit after apply: the 3 (was) lines `% DDF.C1=...`, `% DDF.D11=...`,
   `% DDF.D12=...` restored with `% MMP, 11/21/2026 (was)`: they are the 2025 lines verbatim (option 2,
   commented in 2025), commented code of the prior release. Net: 0 dropped, 4 kept. Removed: entries of
   11/19 (search, maintainer, the bound "about 1e-3"). check: 3 failures, exactly these 3 lines – tool
   limitation, deliberate; codeident OK, w/crlf, numstat OK.
37. `PIETOOLS_examples/examples_PDE_library_PIETOOLS.m` – condense (new notes DJ 5 of 18, of which 4 are 2025
   entries reworded), TAG `DB, MMP, DJ, 11/21/2026` (commits: DB from 05/20/2026, MMP from 10/28, DJ 11/03;
   prep proposed `DB, MMP, DJ, MP`: MP comes from the reworded 2020 Initial coding entry, MP = MMP, not a
   change since 2025). Note (after the 2025 entries, which stay):
   `% DB, MMP, DJ, 11/21/2026: Changes since PIETOOLS 2025: examples 53-55 (3-D)` /
   `%   and 56-72 (mixed domains) added, n_examples = 72; variants 13.1, 13.2 of` /
   `%   the Euler-Bernoulli beam; case 16 (the 4th order Timoshenko beam of the` /
   `%   entry of 05/24/2021) stops with an error: its BCs give an unstable` /
   `%   system; each case collects st, the recommended settings of its example` /
   `%   file; on 'y' settings = st is assigned (cleared without st) and the` /
   `%   executive runs through PIETOOLS_auto_execute_sop; the forms of the tests` /
   `%   are reset with the flags; the 2-D message names st; table of contents` /
   `%   (test, form, stability notion of each case), help and license revised.`
   (The note states, per developer rule 4, that the 4th order implementation of the 2025 entry DJ
   05/24/2021 is no longer called.)
   Hand edits before apply: line 1 marker `% MMP, 10/28/2026` removed (it marked the removal of a stray "2"
   added and removed after 2025; line 1 equals the 2025 code up to a leading blank); the block markers
   `(start)`/`(end)` became `% DB, MMP, DJ, 11/21/2026: BEGIN/END mixed-domain examples 56-72` and
   `... BEGIN/END references [23]-[34]`; the nested pair `(start, example 72)`/`(end, example 72)` removed
   (inside the 56-72 block, development history). Apply 27 / 58 / 99 (prep 27 / 58 / 106; the 7 marker
   lines above were handled by hand). Removed: 13 header entries of 10/28-11/19 (cleaning/ and
   orphaned_file/ paths, maintainer, suite case S17, certification dates, release pass).
   check: 5 failures, exactly the 5 reworded 2025 entries (`MP, 07/01/2020: Initial coding;`,
   `DJ, 05/24/2021`, `05/31/2021`, `12/29/2021`, `08/25/2022`; 2025 wording `Initial coding MP - 7_1_2020`,
   `... DJ - 05-24-2021` etc., reformatted by DJ in 6e65e806, 11/03/2026) – kept per BRIEF step 3, the tool
   does not match the 2025 date formats; codeident OK, w/crlf, numstat OK.
   Open: the case-16 user message (code, not changed) reads `Example 16 was removed from the stability
   examples on 11/19/2026.` – a development date in a release message.
