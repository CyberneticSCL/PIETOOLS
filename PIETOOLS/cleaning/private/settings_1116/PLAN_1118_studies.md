# Settings studies A–F: testing plan and degree-tweaking plan (Sister session, 11/18/2026)

Status: PLAN for approval by the maintainer. Nothing has been run for it. Inputs: HANDOFF_1118.md,
psatz_faces_results/DECISIONS.md (D84–D110), example_settings/EXAMPLE_SETTINGS.md, the library pass
(testing_codes/private/baselines/2026-11-17_tiers/private/results_all_1118.csv), slack_knobs_1d_light_1118.md.
Every number below comes from those records; D-numbers cite DECISIONS.md.

## 0. Rules of every run

- **Tree.** ndopvar at eec4f07d plus the uncommitted working tree of HANDOFF §1. Each run records the
  HEAD hash and the hash of every settings file it reads.
- **Harness.** The 11/17 library harness (lb_run, path 'sop', catalog of 162 rows, size gate
  gate/lpisolve.m at 2e6 variables and m 45,000), copied to a new folder
  testing_codes/private/baselines/2026-11-18_studies/ with one addition: a settings hook (appdata
  `LB_ST_HOOK`, a function of (st, PIE) applied after lpisettings_sop in lb_exec, and an optional
  stability-form override per row). A candidate tier then runs without editing a tier file; with the
  appdata unset the copy runs as committed. The committed 11/17 harness stays unchanged.
- **Records.** One results file per run (results_<study>_<tag>.csv), never results_A/B.csv of the pass;
  summaries under cleaning/private/settings_1116/<study>/.
- **Gates.** Both, per row: harness (|fr − 1| <= 0.3, no numerr, pinf, dinf) and executive (status
  optimal or inaccurate, with eta). A verdict change counts only when both gates agree; rows on which the
  gates disagree are listed apart (D93: 13 light rows of the pass). References: results_all_1118.csv
  (light, heavy) and head_33d30877/ (HEAD 33d30877), compared with lb_compare.pl.
- **Metrics per row.** Verdict; threshold or bound (bisection rows); decision variables and m of the
  certifying pass; number of solves; seconds.
- **Resources.**
  - At most two MATLAB processes on the machine, through the lock files in the matlab_slots folder of
    the Fable scratchpad (one lock per process, sister_<run>.lock).
  - MATLAB sessions under 20 minutes while the maintainer is online, so rows run in batches; longer
    batches only after 04:30.
  - 1-D before 2-D; size-only sweeps before solves.
  - Sol only on the word of the maintainer for any job above ~10 CHE.
- **Unit costs** (one tier over a row group, light, from the 11/17 pass):

  | row group | seconds |
  |---|---|
  | 1-D rows (PDE 1–32, DDE, NDSDDF; 116 rows) | 521 |
  | 2-D library rows (PDE 33–52; 22 rows) | 1,363 |
  | mixed rows (PDE 56–72; 16 rows) | 701 |
  | the same 1-D rows at veryheavy | 11,372 |

## 1. Testing plan by study

### A. D2, the 3-D slack cap (the maintainer: continue)

**Question.** Does some degree choice bring D2 (library case 56 at 0.5 tau*) under m = 45,000 rows and
certify it? If none does, which Sol rung settles it?

**Evidence.**
- The capped light 'Q' form certifies only at no delay (D82).
- Cap 0 with dw_space [1 0] is infeasible with a certificate (D83).
- Memory follows m: MOSEK holds a dense Schur complement at 23–25 bytes per m² (D90). So m = 45,000
  needs about 47–51 GB, and the cap-1 rung (m 82,551) about 160–170 GB.

**Runs.**

| step | content | cost |
|---|---|---|
| A1 | Input pending: d2_size2 (Fable session, size only; not started at 11/18 04:05, waiting for a MATLAB slot; results to cleaning/private/settings_1116/d2/): rungs (cap 1, cap_inc none / 'auto'), (0, 'auto'), (2, 'auto'), (NaN, 'auto'). | none here |
| A2 | Size-only sweep from the smallest rung of A1, one knob at a time: dD ∈ {1, 0}; deg_Q{3} ∈ {1, 0}; psatz{3} ∈ {'faces1', 'auto'}; dw_space ∈ {[0 0], [1 0]}; cap_inc ∈ {[], 'auto'}. Records m, variables and nnz(A) per rung. dD 0 cut the 2-D programs of D99 by 3.1–129x. | 30–45 s per rung; one session of about 15 min |
| A3 | Solve every rung with m <= 45,000 that is not yet refuted, at 0.5 tau*, verdict by both gates; control at 1.1 tau* for any certificate. On the workstation such a solve needs ~50 GB, so it runs alone on the machine (after the Fable agent ends). On Sol: about 15–40 CHE per rung (job 1858048: m 42,222, 30.7 GiB, 498 s of MOSEK, 14 CHE). | maintainer word per Sol job |
| A4 | If no rung under 45,000 certifies: the cap-1 rung (2,721,692 variables, m 82,551; ~170 GB; 7–8 h on 64 cores; ~900–1,000 CHE; D90), or the uncapped program (m 163,847; a highmem node; tens of hours). | maintainer decision; nothing submitted without it |

**Stop rule.** A rung infeasible with a certificate is not raised further in the same knob.

### B. Degree tweaking of the tiers

See Section 2.

### C. The extreme tier (low priority, last)

**Question.** What is the extreme tier for? Today it certifies 41 of 158 rows: 41 infeasible, 42 unknown,
22 lo_fail (D104).

**Steps.**
- **C1 (no runs).** Define the target: the smallest programs that certify the rows light certifies in one
  pass at the library parameters. Those rows are listed from results_all_1118.csv (light, n_solves = 1).
- **C2 (1-D first).** Three candidates on the 116 1-D rows (the extreme 1-D rows cost 78 s today):
  1. as is;
  2. the stripped 1-D entries (storage row, deg_Q{1} [2 2 2]);
  3. candidate 2 with dD 0.

  The best 1-D candidate then runs on the 2-D rows (extreme 2-D: 935 s).

**Deliverable.** Certified count and summed decision variables against stripped and light.

### D. The 'P' executive option in the 1-D delay examples (D100)

**Design (no runs).** The smallest mechanism by which a library block names the stability form:
- **D-a.** A workspace variable `stability_form`, set by the library block and read by
  PIETOOLS_auto_execute_sop in its two stability lines. When absent, behaviour is as today ('Q', 'Qd').
  This edits the auto-executor, which D100 left unchanged.
- **D-b.** PIETOOLS_DDE.m and PIETOOLS_DDF.m call `PIETOOLS_stability_sop(PIE,settings,stability_form)`
  themselves when the block sets the variable, and clear the stability toggle so the auto-executor skips
  it. The auto-executor is not touched; the same three lines live in two entry scripts.
- **Recommendation.** D-a: one place for every entry script, default path unchanged. D-b if the
  auto-executor stays as committed. The maintainer chooses.

**Measurement (through the hook; independent of the design).**
- Rows: the 66 delay rows (41 DDE, 25 NDSDDF) at light in form 'P' ('Pd' on the dual rows), against the
  form-'Q' st blocks committed in eec4f07d.
- Cost: about 8 min per configuration, two configurations, one session.
- Basis: on DDE03–DDE08, light 'P' reaches the form-'Q' thresholds with 1.4–4.6x fewer variables, and on
  DDE03 the documented 1.71785 (D99).
- A row keeps 'Q' wherever 'P' fails to certify what 'Q' with the block st certifies.

**Deliverable.** Per row, the form to name. On approval: the patch (D-a or D-b plus the block lines),
verified by running PIETOOLS_DDE on each changed block.

### E. Library pass gaps

**Rows.**
- veryheavy PDE 56–72 (17 rows)
- stripped PDE 61–72 (12 rows)
- PDE 53–56 at light and heavy (rerun after the lb_build_system guard of D89)

**Procedure.** Clear the keys of those rows from results_A_started.txt and results_B_started.txt first:
53–56 already carry two starts, and a third start records 'crash'. Then run each group into its own file
(results_E_<tier>.csv), size gate on.

**Cost.**
- stripped 61–72 and 53–56 at light and heavy: 15–25 min, two sessions.
- veryheavy 56–72: unpriced on these rows; heavy took 1,328 s on 16 mixed rows. Batches of 3–4 rows
  (<= 20 min each); the size gate refuses the largest programs.

**Deliverable.** results_all_<date>.csv with the rows filled, and the summary table per tier.

### F. The report

**Fragments** (report/):
- what changes: tiertable, knobs, settings_ref, open;
- settings to document: deg_Q, max_raise_Q, raise_Q, wcross, 'faces1', Dmax, cap_inc,
  settings_2d_stock, light deg_Q{1} [4 2 2];
- results and decisions to cite: D91–D108; decisions D100, D109.

**Writing.** Numbers from the records only; writing rules of HANDOFF §3F. Rebuild with MiKTeX (pdflatex
twice, or latexmk), then check the log for undefined references and overfull boxes. No MATLAB.

**Ownership (settled 11/18).** The Fable session confirmed that no report agent of its own was launched
(D109 recorded an intention the handoff superseded); study F belongs to this session alone.

## 2. Degree-tweaking plan (study B)

**Adoption rule.** A tier change is adopted only when it:
1. keeps every certificate of the current tier, by both gates, on the rows of its dimension; and
2. lowers the summed variables or solves, or certifies rows the current tier does not.

Each knob is measured alone before any combination.

### B1. Two-variable lift raise and Q degree (strongest evidence)

**Evidence** (D99, EXAMPLE_SETTINGS §3):
- With dD 0 and deg_Q{2} 0, PDE 33/35/40 certify at 409/8,009/8,518 variables, against
  52,660/72,100/73,463 at light.
- PDE 37 needs deg_Q{2} [1 0 0] (16,656). PDE 34 is infeasible at deg_Q{2} 0.
- PDE 43, 58 and 72 need deg_Q{2} 2.
- dD 0 alone certifies PDE 57 at 140,457 (light 433,725) and PDE 66 at 16,740 (light 56,548).

**Candidates** (light, 2-D entries only, through the hook):

| tag | setting |
|---|---|
| B1.0 | light as is (dD 1, deg_Q{2} 1, raise_Q 'slack', max_raise_Q(2) 1): the control, rerun on the same tree |
| B1.1 | dD 0 in 2-D |
| B1.2 | dD 0, deg_Q{2} 0, raise_Q 'deg', max_raise_Q(2) 2 (Q degree 0 → 1 → 2) |
| B1.3 | dD 1, deg_Q{2} 0, raise_Q 'deg', max_raise_Q(2) 2 |

**Rows.** The 22 2-D library rows and the 16 mixed rows, at light.

**Cost.** About 2,064 s per candidate at the sizes of the 11/17 pass, so about 2.3 h for four. Batches of <= 20 min in
one slot.

**Follow-up.** One run of the winner at heavy (heavy has deg_Q{2} 2).

**Tightness (added 11/18 on the word of the maintainer).** D99 measured verdicts at the library default
parameters, not how tight the bounds are. A candidate is therefore also measured where the answer is
known, and the per-dimension dD decision waits on this data.
- **Stability reach.** For each candidate and for light (dD 1), bisect the parameter toward its known
  margin on seven library cases (bracket [0.5, 1.1] × margin, tolerance 1e-3 of the margin; the
  harness bisection, so 1.1 is a soundness control). Each reach is the largest certified value.
  | case | system | parameter, margin |
  |---|---|---|
  | 34 | 2-D KISS | lam < 2 pi² |
  | 37 | 2-D parabolic | a < 0.5 pi² |
  | 43 | S18 | tau < 1.034726 |
  | 57 | FO2s | tau < 2.881401 |
  | 64 | TL | lam < 13.686938 |
  | 66 | AGE | b < 2.302685 |
  | 72 | S17 | r2 < pi² |

  Implemented as st_run tags T0–T3_<system> (catalog_tight.mat).
- **Gains against exact values.** PDE 45, 46, 47 (exact L2 gains 0.11101, 0.17110, 0.20668), S24 (0.5)
  and fin2 (gap to the record). PDE 45–47 are UNKNOWN in every tier of the pass, on HEAD and on 10/28
  (MOSEK numerr 2; results_all_1118.csv, head_33d30877/). The failure is diagnosed first (form 'Q'
  against 'P' of PIETOOLS_Hinf_gain_sop on PDE 46 at light, one solve each, about 1 min) before any
  candidate is measured on gains.
- **Report per system.** Reach or gain; variables and m of the last certified program; solves; seconds.
  A candidate that cuts the variables but lowers a reach or loosens a gain beyond the solver variation
  (factor 1.3 on gains, D51) is a cost/accuracy trade for the maintainer, not an adoption.
- **Cost.** About 5 min per system and candidate: 8 systems × 4 candidates ≈ 2.5 h in batches of
  <= 20 min.

**Code consequence if adopted.** dD is one number for every dimension in lpi.ineq. A dD for 2-D alone
needs lpi.ineq.dD per dimension (1 x 3), read in lpi_ineq_sop: a small change, made only after the
measurement.

### B2. One-variable lift raise (1-D before 2-D)

**Candidate.** Light with dD 0 on the 116 1-D rows. Cost: about 9 min.

**Adoption.** Only with no certificate lost. In 1-D the lift sets the Markov–Lukacs degree, and the slack
records of D102 and D106 show losses when the slack bases shrink.

### B3. Heavy one-variable Q degree

**Basis.** Heavy and light build the same 1-D 'Q' programs since light moved to [4 2 2] (D94, D96).

**Candidates.** heavy deg_Q{1} ∈ {[5 2 2], [6 2 2]}.
- [5 2 2] reaches the documented DDE03 threshold 1.71785 at 3,337 variables; [4 2 2] reaches 1.7161724
  at 1,712.
- [6 2 2] reaches it at 4,174.
- Both reproduce the [4 2 2] thresholds on DDE04, 05, 07, 08 (EXAMPLE_SETTINGS §5.1).

**Rows.** The 66 delay rows and PDE 1–32, at heavy. Cost: about 8 min per candidate.

### B4. Veryheavy one-variable Q degree

**Candidate.** deg_Q{1} [8 6 6]. On DDE08 it reaches 3.9144 at 30,356 variables, against 3.8957 at 50,214
for [6 6 6]; the documented range is 3.9–3.95.

**Rows and cost.** Veryheavy 1-D costs 3.2 h per pass (DDE17 alone about 50 min, D94), so only the delay
rows run, after 04:30.

### B5. Triangular cap per example (no runs)

From slack_knobs_1d_light_1118.md, list the 1-D rows where cap_inc 'auto' lowers the variables at the
same verdict and threshold. These are candidates for the st blocks of examples (D106), not for the tiers.

### B6. Three variables

Study A.

### Not proposed (already decided)

- 1-D degrees of light, stripped and custom (D109).
- max_raise_Q(1) 0 as a tier change: example-specific (D109).

## 3. Order and budget

1. **Now** (one slot while the Fable agent holds the other):
   - A2, once d2_size2 is posted (about 15 min);
   - B2 (9 min);
   - B3 (2 × 8 min);
   - the D measurement (16 min);
   - E, stripped and 53–56 (15–25 min).
2. **After the Fable agent ends** (about 06:10): B1 in batches (about 2.3 h). A3 solves near m 45,000
   need the whole machine, or Sol.
3. **After 04:30 of the next night:** B4 and E veryheavy.
4. **F:** in parallel; no MATLAB.
5. **C:** last.

## 4. Decisions asked of the maintainer

1. The settings hook in a copy of the harness (test code only).
2. Gate for verdicts: both gates, a change counted only when both agree (proposed), or the executive
   eta gate alone (D93).
3. D: mechanism D-a or D-b.
4. A: each Sol rung, priced before submission.
5. If B1.1 or B1.2 wins: a per-dimension dD (code change in lpi_ineq_sop and the tier files). Maintainer
   (11/18): first more data on tightness, the reach and gain measurements added under B1.
6. Whether D109 (max_raise_Q 0 is example-specific) also covers veryheavy, where the ladder spends four
   solves per infeasible step (DDE17 about 50 min, D94).

## Appendix: for the eye of the maintainer (HANDOFF §3G, no study)

- The DEMO12 witness reports a closed-loop gain of 5.6486 at omega 7.8e5, above both certified bounds:
  synthesis 5.6253 and analysis 5.6237 (analysis eta 1.5e-5).
- Four tests under testing_codes/private call poslpivar_settings_sop with tier settings and would fail.
