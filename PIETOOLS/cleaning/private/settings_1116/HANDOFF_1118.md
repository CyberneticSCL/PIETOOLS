# Handoff of the settings studies to the Sister session (MMP, 11/18/2026 04:15)

Maintainer, 04:12: "Let move the outstanding studies to Sister agent and save you for Fable-level
issues which may arise later. Provide enough context to the sister agent to propose a testing and
degree-tweaking plan." The first deliverable asked of the Sister session is therefore a PLAN (a
testing plan and a degree-tweaking plan) for the maintainer's approval, not runs.

## 1. State of the repository (branch ndopvar, all pushed)
- bb367fe1: the settings revision on the layout of the report of 11/16 (356 files). Tiers on a
  storage row per dimension; lpi.ineq (dD, dw, wmin, wcross 0, psatz {'auto','faces1','auto'},
  Dmax false, cap_inc []); lpi.var.deg (synthesis free operators only); sop (slack 'like',
  dw_feas, dw_P, dw_Q, max_raise_feas, max_raise, cap_delta, margin, and deg_Q, max_raise_Q,
  raise_Q for the direct stability form 'Q'). settings_2d left the tiers for
  settings_2d_stock(level). Light's 1-D deg_Q is [4 2 2]; stripped and custom keep [2 2 2]
  (maintainer: fine). lpisolve resolves the tiers' empty solver field.
- 03a8a261: DEMO12 (Dan Braghini's inlined synthesis) runs on the tier layout.
- eec4f07d: example-specific settings in fourteen examples (the Opus agent's st blocks; its note
  psatz_faces_results/example_settings/EXAMPLE_SETTINGS.md).
- Uncommitted: get_lift_degs with the lift ceiling moved after the weight rule (knob off); the
  records (psatz_faces_results/DECISIONS.md D84-D109, mixed_domain/GAP_RUN_LOG.md), the notes in
  this folder (newton_polytope_guidance.md, q_p_differentiation_proposal.md,
  slack_knobs_1d_light_1118*.md, DEMO12_report_1118.md, commit_plan_1118.md), the final rows of
  the library pass (testing_codes/private/baselines/2026-11-17_tiers/private/results_B.csv,
  results_all_1118.csv, final_summary_1118.md). Commit them with the next batch on the
  maintainer's word (he gives it per commit and per push).

## 2. What is established (cite the DECISIONS entry)
- D91-D94: the library pass under all six tiers against the 10/28 reference and a HEAD reference
  run (head_33d30877/). The commits between 10/28 and HEAD change no DDE/NDSDDF value; the
  working tree lowered five light DDE thresholds by 0.3-5.3% because the direct form reads Q from
  the light entry [2 2 2]; [4 2 2] (now set) restores heavy's thresholds on all ten 1-D bisection
  rows (results_ddeQ4.csv). In 2-D the direct form certifies seven rows the storage form did not
  and cuts PDE 33-37 heavy to a third of the variables; PDE 39 light is MOSEK OPTIMAL (eta 2.4e-9)
  but OPT_STATUS 7.14, refused by the harness gate (|fr-1| <= 0.3) and accepted by the executive's
  eta gate; the term set is not the cause; Q at the heavy 2-D degree costs 3.7M variables.
- D95-D96: the separation of the direct form (deg_Q, max_raise_Q, raise_Q); verified bit-identical
  with the fallback (ab_sep.m). The 'P' form at light certifies DDE05 at 1.5707 with 651 variables
  against 1,549 for 'Q' at [4 2 2].
- D99: the example-specific settings (14 examples): PDE 33/35/40 at 409/8,009/8,518 variables
  against light 52,660/72,100/73,463 (dD 0, deg_Q{2} 0), at the library default parameters only;
  DDE08 3.9144 at 30,356 variables with Q [8 6 6]; in 1-D the form 'P' at light beats every 'Q'
  setting by 1.4-4.6x in variables; in 2-D 'Q' wins ('P' infeasible on 43, 57, 58). Minimal Q
  degrees reaching light-P's verdict: DDE03 [5 2 2], DDE04 [3 2 2], DDE05/07/08 [4 2 2]; 2-D 0 for
  33/35/40, [1 0 0] for 37. deg T reads 0 on the DDE rows, so the containment rule deg Q >= deg P +
  deg T does NOT explain the light loss (D92's "T of degree 1" withdrawn).
- D100: NO form field in the settings: 'P' and 'Q' prove different forms of stability (PDE-level);
  the example file names the executive option; 'P' is used wherever it certifies (the stronger
  form). The auto-executor is NOT changed (maintainer, 03:00).
- D102, D106: Shivakumar's bounds as knobs (newton_polytope_guidance.md): the triangular cap
  (cap_inc 'auto') saves 17% of the variables on the 124 1-D rows at light and loses 8 of 124
  certificates (DDE01 dual, DDE13, DDE18, PDE 21 dual, 30, 31, 51, 52); the lift ceiling with the
  weight held costs 2.42x for bounds within 0.1% (4 numerical losses). Both knobs off in the tiers;
  per-example candidates. The derivative targets are purely integral on L2 (maintainer), so the
  bounds apply to the slacks that set the cost.
- D109: max_raise_Q 0 in 1-D is an example-specific special case, not a tier change; the extreme
  tier revision is low priority and needs a study; the D2 investigation (3-D cap) continues; the
  converter edit is not expected to solve anything (conversion <= a few GB against a 160-600 GB
  factorization: shelved); the report pdf is to be updated with the tier table.

## 3. Outstanding studies (propose a plan; run on approval)
A. D2 (the 2-D delayed heat equation with its 3-D delay state, library case 56, tau = 0.5 tau*):
   capped light 'Q' certifies only at no delay (D82); cap 0 + dw_space [1 0] infeasible with a
   certificate (D83); uncapped needs ~600 GiB (23-25 bytes per m^2 in MOSEK's Schur complement,
   D90); the cap-1 rung is 2,721,692 variables, m 82,551 (~170 GB, 7-8 h on 64 Sol cores, ~1,000
   CHE: Sol jobs above ~10 CHE need the maintainer's word). A size-only run of the cap rungs with
   the triangular cap (cap_inc 'auto') is in progress on the workstation (scratchpad
   sol1117/size/d2_size2.m and .log of this session; copied to d2/ here when it ends): the
   question is whether cap_inc or another degree choice brings a rung under ~45,000 rows.
   Tooling: sol1117/ (scripts, shadow lpisolve that sizes without solving), the Sol folder
   /scratch/mpeet/pietools/gap1117 (memory file asu-sol-access-and-environment.md).
B. Degree tweaking of the tiers (a plan first): candidates with evidence are per-dimension
   deg_Q (light 1-D at [4 2 2] done; 2-D entries 1 at light are below the containment rule but the
   rule is unaffordable there), max_raise_Q (example-specific), cap_inc per example, stripped and
   custom (maintainer: fine as they are). Any tier change must be measured through the harness
   (testing_codes/private/baselines/2026-11-17_tiers/private: lb_run path 'sop', the catalog of
   162 rows, lb_compare.pl against results_all_1118.csv and head_33d30877/) with the size gate
   (gate/lpisolve.m, 2e6 variables, m 45,000) and two MATLAB processes at most (lock files in the
   scratchpad matlab_slots of whichever session runs; sessions under 20 min while the maintainer
   is online). Test 1-D before 2-D.
C. The extreme tier: low priority; a study of what it should be (today: the old degrees on the
   new layout, patterns {'a'}, {'aa'}, {'aaa'}, wmin 0, max_raise 0).
D. The 'P' executive option in the 1-D delay examples (D100): the DDE library blocks set st
   (eec4f07d) but cannot name the form; design the smallest mechanism for an example to name it
   (e.g. a variable the library block sets and PIETOOLS_DDE.m / the auto-executor reads) and
   measure the 1-D delay rows in form 'P' at light against the committed form-'Q' st.
E. Library pass gaps (postponed by the maintainer): veryheavy rows 56-72, stripped 61-72, cases
   53-56 at light and heavy (results_A_started.txt counts two starts for those keys: clear them
   before a third start, else lb_run records 'crash').
F. Documentation: update the report (report/settings_report_2026-11-16.tex with its fragments,
   tier table in settings_report_tiertable.tex, knobs in settings_report_knobs.tex, reference in
   settings_report_settings_ref.tex, open items in settings_report_open.tex) for deg_Q,
   max_raise_Q, raise_Q, wcross, 'faces1', Dmax, cap_inc, settings_2d_stock, light's [4 2 2],
   the measured results D91-D108 and the decisions D100, D109; rebuild the pdf with MiKTeX
   (pdflatex/latexmk on the path). Writing rules: no opinions (every comparative with its metric),
   no "we" that excludes the reader, active voice, no banned words, spaced en dashes, name every
   Positivstellensatz term by its factor (plain term; box product; the faces (s_d - a_d),
   (b_d - s_d); box quadratics; face products), numbers from the records only.
G. For the maintainer's eye, not a study: DEMO12's witness reports a closed-loop gain (5.6486 at
   omega 7.8e5) above both certified bounds (synthesis 5.6253, analysis 5.6237; analysis eta
   1.5e-5 above the harness gate); four tests under testing_codes/private still call
   poslpivar_settings_sop with tier settings and would fail.

## 4. Running now
- The example-settings Opus agent of the Fable session (resumed 04:05, about 2 h): PDE 36, 45-48
  and the size-gated 50, 61, 65, 69, 70; it reports to the Fable session, which relays and commits
  on the maintainer's word. It holds one MATLAB slot at a time (lock files exst_*.lock in the Fable
  scratchpad matlab_slots): coordinate through the same folder, two processes at most on the
  machine.
- d2_size2.m (size only, minutes).

## 5. Rules that bind every session
Annotation protocol (header entry, column-77 markers "MMP, <date>", (was) lines; terse comments);
report other developers' defects, do not patch them (Declan, Dan Braghini); results, logs and
reports under cleaning/private/<topic>/ (removed before release), tests under
testing_codes/private/; SOSTOOLS files never mention PIETOOLS; the maintainer's word for every
commit, push and Sol job above ~10 CHE; one results file per run, never results_A/B.csv of the
pass; verify each direction of an inverse pair against its semantics; measure, do not reason
about cost (time and memory over a sweep).
