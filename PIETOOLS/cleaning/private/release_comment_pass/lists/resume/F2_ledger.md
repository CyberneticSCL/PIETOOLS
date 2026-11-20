# Ledger – F2 (11/20/2026)

Lists: `F2.tsv` (79 files under PIETOOLS_examples/; 57 new = 0, 22 new = 1) and `F2_runs.txt` (7 runs in 5 files;
the PDE library is on both lists). 83 files read and edited, comment-only. Final checks (after HEAD moved to
41c9210c, which touched none of these files):

- `codeident.sh` on all 83: **83 OK, 0 DIFF, 0 NEW**.
- `runs.awk` (N = 20) over all 83: **0 runs**. Runs of 15–19 left by judgment: examples_DDE_library_PIETOOLS.m
  714-728 (15; matrix data under a comment giving the equations), examples_PDE_library_PIETOOLS.m 2031-2048 (18;
  the if/elseif chain of process_inputs under a comment naming the branches).
- Line endings and final-newline state of all 79 F2 files equal the state recorded at the start (77 CRLF, 2 LF:
  2D_Plate_Eq, 2D_Wave_Eq_Damped; git warns that autocrlf would convert those two). The 4 run-only files keep CRLF
  and their HEAD final-newline state.
- All 79: GPL block once, support line once, an Initial coding line or entry; the 57 new = 0 files carry
  `% MMP, 11/19/2026: Release comment pass (comment-only).` once; the 22 new = 1 files have no dated log entry,
  per-line marker or (was) line left.

Process note: the three libraries were edited by a forked agent of this session that stopped at the usage limit
(01:50) before writing its notes; its edits (codeident OK, EOL kept) were reviewed from the diff and NDSDDF was
finished here. Common to every example file below unless stated: a what-it-does paragraph in the
`% [OUT] = NAME(IN) ...` form naming the test flags, the form (stability/gain 'P' = storage <T x, P T x> with
P >= eppos I, 'Q' = generalized storage T'Q, duals 'Pd'/'Qd'; H2 'o'/'oco' and duals) and the GUI behaviour;
license; support line; authors from git (Declan -> DJ, DanBraghini -> DB, Sachin -> SS, MATTHEW M PEET -> MMP);
step comments for the default parameters and the GUI branch; a comment stating what the settings set (tier,
storage lengths/patterns, slack lift degree dD, cap_inc, psatz set, wcross, deg_Q, free-operator degrees,
weight raises).

## Edited files

Libraries (new = 0)
- examples_PDE_library_PIETOOLS.m: what-it-does; INPUT corrected to process_inputs (index range 1-72 with
  variants, numeric format codes 1/2/3, random example for a format alone, outputs, which examples have batch
  and GUI files); pointer to PIETOOLS_PDE.m removed; license; support; release line; a comment on the switch
  explaining the forms; a `Test:` line (flags and form) at every case; ToC entries matched to the example files
  (Das -Bxr and threshold 9.7946; Tang x(0) = -xo; Ahmadi_2 Cm; Timoshenko_1 states and c; Timoshenko_2 c = 0
  stable; 2-D advection c1 = c2 = -1; 2-D ODE term x2_{(2,0)}; reaction-diffusion r; Krstic; x2(1+tau) = u;
  PDDE states x2, x3 over [0,pi]; Unstable control +s*u; controller outputs; Demetriou [16]); tier comparisons
  and solver remarks removed (cases 20-24, controllers, H2 cases, PDDE); case 16 comment in present tense
  without the maintainer quote and orphaned_file path; cases 100-101 stated unreachable; process_inputs header
  (what/INPUT/OUTPUT) and step comments (breaks the 37-line run 1890-1926). Here also: gains in the ToC set to
  the verified values, case 21 "(L2 gain 0.8419 for k = 0.5)", case 23 "(gamma = 23.706 ...)", case 24
  8.2136 -> 8.2141 (ne = 1).
- examples_DDE_library_PIETOOLS.m: what-it-does with the DDE equations; pointer to PIETOOLS_DDE.m removed;
  license; support; release line; a settings/form description before every st block; test numerics removed
  (gam/IPM lines, "minimal comp: 154.95s", "requires veryheavy", "'P' at light: 1.42229"); the [5] block comment
  without history framing; step comments in the active block (Example 2b of [13]).
- examples_NDSDDF_library_PIETOOLS.m: header name corrected (was examples_DDF_library); what-it-does with the NDS
  and DDF equations; pointer to PIETOOLS_DDF.m removed; license; support; release line; settings/form
  description before each of the 14 st blocks; dims/cpusecs and N=5/N=30 gam/IPM lines removed; the 11/19
  comment of the active block in present tense; step comments in the active showering-people block.

2D folder (new = 0 unless stated)
- Hinf_Gain_Examples/..._2D_ReactionDiffusion_DDDD.m, _DDDN.m: ne added to the parameters; settings (heavy;
  2-D storage lengths/pattern; dw_P 1).
- Hinf_Gain_Examples/..._2D_ReactionDiffusion_DNDN.m: gain series corrected to mu_{m,n} (was mu_{2m-1,2n-1});
  the code and the stated gam = 0.20668 use mu_{m,n} (both series evaluated: 0.20668 vs 0.20609); settings.
- Hinf_Optimal_Observer_Examples/..._2D_ReactionDiffusion_DNDN_EdgeObserve.m: outputs corrected to the code
  (y1 = int x ds2 + w2, y2 = int x ds1 + w3; the header had the edge values that only the commented term-based
  block uses, now noted); settings (heavy).
- PIETOOLS_PDDE_Ex.m: description matched to the code (states x2, x3; output over [0,pi]; ranges of s1, s2;
  delays 0.25, 0.3); custom-settings H2 bounds removed; license between description and log (DB - 01/20/2025);
  tau-after-params stated; step comments; settings (light; dD 0).
- ..._2D_Heat_Eq_with_ODE.m, _2D_KISS_Model.m, _2D_NS_Antonelli.m, _2D_Parabolic_Eq.m,
  _2D_Reaction_Diffusion_Eq.m, _2D_Telegraph_Eq.m, _2D_Transport_Eq.m, _2D_Wave_Eq.m,
  _Coupled_ODE_Heat_Eq_w_Delay.m, _Heat_Eq_w_Delayed_Boundary_Input.m, _Heat_Eq_w_Interior_Delay.m,
  _Wave_Eq_w_Boundary_Delay.m, _Wave_Eq_w_Boundary_Delay_2D.m: PDE_b removed from OUTPUT (not returned);
  "(header)" pointers in the st entry replaced; settings comments. Per file: Heat_Eq_with_ODE x_{(2,0)} ->
  x2_{(2,0)}, c1, c2 added; KISS/Telegraph ne added; NS_Antonelli states named, nu added; Reaction_Diffusion
  exponential stability for C > 0 (eigenvalues) added; Telegraph (C^2/2*lam) -> (C^2/(2*lam)); Transport
  defaults and finite-time extinction stated; Wave_Eq energy conservation (stable, not exponentially stable)
  stated; Delayed_Boundary_Input "Kristic" -> "Krstic", open-loop threshold lam < pi^2, expand_delays comment;
  Interior_Delay st entry without "certifies"/settings report (3 marked lines reworded), domain [a,b];
  Wave_Eq_w_Boundary_Delay "s1" -> "s" (header, inline comment); _2D variant a, b with the s1 = 1 caveat.
- ..._2D_Plate_Eq.m (1), ..._2D_Wave_Eq_Damped.m (1): log (CC 09/28 entry with its continuation, MMP 11/19),
  markers, (was) lines removed; Initial coding MMP - 09/28/2026 (git author; the log said "CC"); references to
  a demo file and the retired 2-D stability executive removed; "the default stability test cannot certify this
  system" removed; paper results kept without solver details; "no executive flag" restated.

3D (new = 1)
- ..._3D_Reaction_Diffusion_Eq.m, _3D_Advection_Diffusion_Eq.m, _3D_Coupled_Heat_Eq.m: log, markers, (was)
  lines removed; Initial coding DJ - 11/02/2026 (the files' own entry); PDE_b removed from OUTPUT (Reaction);
  thresholds 14.8044, 30.4838, 16.2174 checked; settings (stripped).

Hinf_Gain_Examples 1-D (new = 0)
- ..._Diffusive_Eq_with_Distributed_Disturbance.m (gamma 0.8102 checked by modal series), _2.m,
  _Transport_Eq_with_Disturbance.m: settings; dual form named where set.
- ..._Heat_Eq_with_Distributed_Disturbance.m: license between the description and the log; settings.
- ..._Parabolic_Eq_with_Disturbance.m: gain corrected "15.147" -> "23.706, attained at frequency 0"
  (finite differences N = 1000-4000, converged to 23.7058; frequency grid maximum at 0).
- ..._Reaction_Diffusion_Eq_with_Disturbance.m: gain corrected "8.1069" -> "8.2141" for ne = 1, with the closed
  form (2 tan(sqrt(lam)/2)/sqrt(lam) - 1)/lam (A self-adjoint, C = B*: peak at frequency 0); ne added.
- ..._Wave_Eq_Tip_Damped.m: gain corrected "2" -> "0.8419" for k = 0.5 (|G(i w)| of the transfer function,
  w <= 40, peak 0.84188 at w = 1.665).

Controllers, observers, Examples_Library root (new = 0)
- ..._Euler_Bernoulli_Beam_with_Control.m: output corrected to the code, z = int (1-s) v2 ds.
- ..._Reaction_Diffusion_Boundary_ODE_Control.m: settings.
- ..._Reaction_Diffusion_Eq_with_Controlled_Input.m: titles "pure transport equation" and "with z=0 and w=0"
  replaced by "Reaction-diffusion equation with distributed control u".
- ..._Wave_Eq_Boundary_ODE_Control.m: output order corrected to the code, z = [int x1; xo].
- ..._Heat_Eq_with_Observed_Output.m: settings. ..._Transport_Eq_with_Observed_Output.m: PDE corrected to the
  code, x_t = x_s + w (header had (s-s^2) w).
- ..._Transport_Eq_with_Disturbance2.m: PDE_t entry rejoined (split by the st lines); license between
  description and log (DB - 01/06/2025).
- ..._Unstable_ReactionDiffusion.m: PDE_b removed; "Parameter c can be set" -> "No parameter can be set".
  ..._Unstable_ReactionDiffusion_Control.m: control term s*u added to the header PDE; PARAMS unused and the GUI
  error stated.

Mixed_Domain_Examples (new = 1, 17 files): log entries, markers and (was) lines removed (Initial coding
MMP - 11/17/2026 from each file's own entry); st entries without "certifies"/settings-report references;
provenance tags [read]/[derived]/[abstract] kept. Per file: Kang (light; dD 0; deg_Q [0 1 1]); 2D_Delayed_Heat_Eq
_Fridman "library example 43" removed; Lidstrom finite-difference check values removed; Damped_Wave_Interior_
Delay_Fridman and _Riemann_ cross-file pointers removed; Nicaise, Heat_Edge_Cascade, Road_Field, Two_Layer
thresholds/abscissae kept at the finer grid; Heat_Laplace_Beltrami finite-difference check removed;
Transport_Edge_Diffusion "mixed-domain test suite, system S17" and the unchecked claim that no coercive P
certifies it removed; Wave_Delayed_Boundary_Damping_Xu "Library example 44" removed; settings comments on all.

Stability_Examples (new = 0, 19 files) and command_line_input_Examples/stability_15.m (new = 0)
- All 19: settings comments; dual forms named. Diagne: two step comments in the parameter block. Euler_Bernoulli
  _Beam_Eq: INDEX added to INPUT (13/13.1 vs 13.2), energy conservation stated. Das: the 11/19 header note
  without its run date. Tang: BC corrected to the code, x(0) = -xo. Lamare, Saba: INDEX added to INPUT.
  Parabolic_Eq_Gahlawat: parameter names corrected (a, b, c, not afun, bfun, cfun) and c computed before params
  stated. Saba: "% Stable with stripped settings" removed from two parameter lines. Timoshenko_Beam_1: damping c
  added to the parameters. Transport_Eq: v and finite-time extinction stated. Datko_2: parameters stated.
- stability_15.m: STABILITY_15 description (system as coded: no damping, c = 0; displayed; workspace PDE, PIE;
  params); "(unstable)" corrected twice to "c = 0: stable, not exponentially stable" (the 11/19 log entry of the
  file states this after the x4_{ss} restoration); license; authors SS - 11/17/2022; release line.

Runs (`F2_runs.txt`; 4 files passed on 11/19, codeident OK before editing)
- snippets_from_manual/PIETOOLS_Code_Illustrations_Ch2_Introduction.m (126-155, 209-238): panel comments.
- .../Ch4_Declaring_PDEs.m (74-97; also 47-64 at 17): step comments.
- .../Ch6_PIESIM.m (150-195, 240-268; also 96-110 at 15): panel/label comments.
- command_line_input_Examples/stability_3.m (64-83): two step comments in the parameter block.
- examples_PDE_library_PIETOOLS.m (1890-1926): see the libraries.

## Skipped
None (all 79 F2 files clean in git status before their first edit; the 4 run files codeident OK).

## Kept log entries citing cleaning/, testing_codes/, test numerics, maintainer quotes or run dates (new = 0)
- examples_PDE_library_PIETOOLS.m: MMP 11/08 "(cleaning/private)"; MMP 11/17 "S17 of the mixed-domain suite",
  "settings that certify the example"; MMP 11/18 "measured against the tiers in cleaning/private/settings_1116/
  psatz_faces_results/example_settings"; MMP 11/19 "(maintainer: no default settings)", "(certified
  11/19/2026)", case 16 "(maintainer)", case 8 "(maintainer: a typo in the paper)".
- examples_DDE_library_PIETOOLS.m: MMP 11/18 (two entries) cleaning/private/settings_1116/.../EXAMPLE_SETTINGS.md,
  taumax values, "1.65 times the variables", "11 SDP solves (heavy: 15-29)", "1.4-4.6 times fewer variables",
  D100, "'P' at light 1.42229"; MMP 11/19 "per-example search of 11/19/2026 (maintainer: no default settings)".
- examples_NDSDDF_library_PIETOOLS.m: MMP 11/19 "per-example search ... (maintainer: no default settings)";
  MMP 11/19 "the bound fell to about 1e-3".
- 53 example files (every new = 0 example file except stability_15): MMP 11/19 "from the (per-example) search of
  11/19/2026 (maintainer: no default settings)" or "no measured setting certifies it; (maintainer: every example
  names its settings)".
- 2D_Heat_Eq_with_ODE, 2D_Parabolic, 2D_Reaction_Diffusion, 2D_Transport: MMP 11/18 "certify the library default at
  fewer decision variables than the light tier"; 2D_Wave_Eq: "... which the light tier does not, at fewer
  decision variables than heavy"; Heat_Eq_w_Interior_Delay: MMP 11/17 "the light tier certifies it in neither
  form", MMP 11/18 "half the variables".
- Heat_Eq_with_Distributed_Disturbance: MMP 11/05 "the H-infinity LPI returns 0.20847".
- Heat_Eq_with_ODE_Das: MMP 11/19 maintainer quote and "Chebyshev, 56 and 72 nodes".
- Saba: MMP 11/19 "maintainer: 'Apply all four'", "real eigenvalue 1.682". Timoshenko_Beam_2: MMP 11/19 commit
  72de7c6b, maintainer quote. stability_15: MMP 11/19 cites PIETOOLS_PDE_Ex_Timoshenko_Beam_2.m.

## Literature citations kept
Bibtex/references in the files (Holmes 1994; Antonelli 2021; Demetriou 2019; Kang & Fridman 2017; Krstic 2009;
Caliskan 2009; Xu et al. 2006; Valmorbida 2014; Shivakumar 2019, 2022; Peet 2019, 2021; Das 2018; Diagne 2012;
Lamare 2016; Saba 2019; Tang 2011; Gahlawat 2017; Ahmadi 2014, 2015; Datko 1986; Jagt & Peet 2025; Fridman &
Orlov 2009; Lidstrom et al. 2016; Turing 1952/Murray 2003; Beale 1976; Kang 2022; Marzocchi et al. 2002; Nicaise
et al. 2009; Berestycki et al. 2013; Hickson et al. 2009; Batty et al. 2016, 2019; Goldstein 2006; Vazquez &
Vitillaro 2011; Zhang & Zuazua 2004; Mazanti et al. arXiv 2002.06128; NDS/DDE library reference lists).

## Candidates for the maintainer (not fixed)
1. ..._Euler_Bernoulli_Beam_Eq.m and ..._Saba.m: the error branch uses `indx`, undefined (the intended message
   never prints; MATLAB reports an undefined variable).
2. ..._Euler_Bernoulli_Beam_with_Control.m: z = int (1-s) v2 ds has no control penalty (D12 = 0); the header
   previously stated z = [u; int 0.5(1-s)^2 v2 ds]. An H-infinity synthesis without a u-penalty is degenerate.
3. ..._Unstable_ReactionDiffusion.m: GUI = 1 loads PIETOOLS_PDE_Ex_Euler_Bernoulli_Beam_Eq_GUI.mat.
4. ..._Unstable_ReactionDiffusion(_Control).m: no params loop (PARAMS ignored); both set st.sop.deg_Q, a setting
   of the stability form 'Q', in H2 synthesis examples.
5. ..._Parabolic_Eq_with_Disturbance.m, ..._Reaction_Diffusion_Eq_with_Disturbance.m print "With sedumi solver,
   this example requires at least settings veryheavy" while st is stripped / light.
6. ..._Parabolic_Eq_Gahlawat.m: c is computed from lam before params, so a params override of lam leaves c
   unchanged. ..._Diagne.m: the first params loop runs before the defaults and is overwritten.
7. PIETOOLS_PDDE_Ex.m: tau is set after the params loop (params cannot change it); pvar s t unused.
8. ..._Wave_Eq_w_Boundary_Delay_2D.m: a, b are parameters, but the BCs evaluate x2, x3 at s1 = 1.
9. ..._2D_ReactionDiffusion_DNDN_EdgeObserve.m: file name and commented term-based block observe the edges; the
   active code senses integrals over s2 and s1.
10. DDDD, DDDN, DNDN: r_lim and L2_gain are computed and discarded (not returned or displayed).
11. stability_15.m: `clc; clear;` runs before `if exist('params','var')`, so params never applies (as listed for
    stability_2-14 on 11/19); the defaults k, aa, II, g are unused.
12. Thirteen 2-D files listed PDE_b in OUTPUT without returning it (comments fixed only).
13. ..._Transport_Eq_with_Disturbance2.m: pvar s1 s2 unused. Ch6 snippet: x3 plot title shows s_{3} for a
    function of s2 (string, not changed).

## Unresolved
1. ..._Diffusive_Eq_with_Distributed_Disturbance_2.m states gamma = 4.23 for R = 21-1e-3; the PDE library ToC
   (case 26) states 2.1457. Not checked (3-state system with s-dependent Cm).
2. Estimator values differ between files and ToC: Transport_Eq_with_Observed_Output "Answer: 1.0012" vs ToC
   "(gain: 1.0004)"; Heat_Eq_with_Observed_Output "Answer: 1.0045" vs "(gain: 1.0005)". Left as written.
3. Ahmadi examples: "Stable for R<2.7" / "R<21" with defaults R = 2.9 / 30 (the ToC calls the bounds not tight);
   left as written.
4. The fork's ToC corrections in the PDE library were checked against the example files edited here; its header
   statements on batch/GUI availability per example (1-15, 17-29; GUI files for 1-15, 17-28, 51) were not
   re-derived.
