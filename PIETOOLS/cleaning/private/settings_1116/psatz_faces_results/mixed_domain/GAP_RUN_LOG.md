# Autonomous gap run, 11/17/2026 (12 hours from 04:35)

Maintainer instruction (11/17, 04:30, offline afterwards): "Run autonomously for 12 hours to find
ways to eliminate remaining gaps, including variations of psatz terms and degrees to match the
structure of the operator inequalities. Explore both Q and P forms. Use subagent for Opus-class
tasks." Earlier the same hour: application-specific settings files, or modifications of the
settings in the example files, are allowed where the default tiers do not suit an application;
settings_2d leaves the proposed structure (its only reader is the 2-D branch of exec_tools_sop,
which work package 4 routes through the row).

## The gaps (stable systems the tiers do not certify; state at 04:30)

| id | system | point | storage form 'P' | direct form 'Q' |
|---|---|---|---|---|
| G1 | S17 transport-edge diffusion | 0.9 r2* | uncertifiable, every storage (proven 11/16) | 0.5: certified (faces, scalar offset); 0.9: not with faces at dw 0 |
| G2 | S18 interior delay (Fridman-Orlov 2009 Ex. 1) | 0.9 tau* | light, heavy fail; veryheavy 2-D entry (Rlight2) certifies at 1,512,873 vars (report tab:bl-mixed) | closed at light with dw 1 and the faces at the scalar offset, 440,023 vars (tab:t4) |
| G3 | age-structured population (lit. item 16) | 0.9 b* | light infeasible 187,984; heavy certifies 751,849 | light certifies: facesS 56,548, faces0 210,340 |
| G4 | delayed damped wave (lit. item 9) | FO2r any tau; FO2s > 0.916 tau* | FO2r infeasible even at 0.1 tau*; FO2s infeasible at 0.5 | FO2s faces0 reach 0.916, numerr above; FO2r infeasible |
| G5 | 3-D heat | 0.9 kappa* | light infeasible, every storage | faces0 2.71M vars, 759 s (C-series) |
| G6 | 2-D delayed heat, 3-D registry (lit. item 8) | 0.5 tau* | infeasible 1.67M, no raise | faces0 47.4M vars, refused |

Not gaps: S3, S8 (the row certifies in 'P' where stock fails; WP4 after-run), items 1, 2, 3
(not exponentially stable), items 5, 10, 11, 12, 13, 15 (certified), item 6 (closed to 3.9e-4).

## Program

- Phase A (Opus agent, S/gapstudy): the knob ladder in both forms on G1 to G6: term sets (none,
  facesS = 'auto', faces0, faces at offset 1, box, corners, product2, mixedfb, offsets 0/1/2),
  dw 0/1/2, wcross 0, Q degrees 1 and 2, storage rows and tiers, cap_delta, max_raise_feas.
  Record: scratchpad gapstudy/GAP_RESULTS.md, copied here as gap_results.md.
- Phase B (design here, port and runs by an Opus agent, S/gapstudy_b): structure-matched slack
  degrees: get_lift_degs per space (D_{k,d}, w_{k,d} from the blocks touching space k, not the
  registry maximum), per-space dw (dw_space), the 'auto' offset switch on the least w over
  spaces and directions. Record: perspace_results.md.
- Phase C: closing configurations into the example files (Examples_Library/Mixed_Domain_Examples)
  as application-specific settings, or a tier change when it loses nothing on the baselines;
  rerun the baselines of the report under any tier change.
- Phase D: report (application-specific settings subsection, open items, gaps closed), DECISIONS
  D66 onward, memory.

Constraints: at most two MATLAB -batch processes (lock files in scratchpad matlab_slots/); the
size cap 3,000,000 variables (litbuild lpisolve shadow); watchdog 1,800 s per script; no
repository code edits by the agents; commits only on the maintainer's word.

## Log

- 04:35 Phase A launched (gap-study agent). WP4 after-run in progress (its agent resumed to write
  the before/after table). Report agent resumed: macro fix, settings_2d out of Section 5,
  subsection "Application-specific settings". Examples-library agent running. Litbuild item 4
  rerun running.
- 05:00 Examples library done: cases 56-71 in Examples_Library/Mixed_Domain_Examples, index raised
  to n_examples 71 (uncommitted); three executive checks match the litbuild record.
- 05:06 WP4 after-run cut by hand on its last rows (S17 'Q' stock at one raise, 5.24M variables);
  WP4 table in scratchpad wp4/WP4_RESULTS.md: the direct 'Q' form certifies S3, S8, S18(0.5) at
  1.19M, 2.32M, 1.33M (3-6x the harness) and fails S17 at 0.5 (UNKNOWN 1.31M): sized "like P".
- 05:25 Executive fix (MMP): the 'Q' form sizes both slacks from their own support (no 'like'),
  as the harness; header entry updated. WP4b agent reruns the 'Q' rows with dw_feas as is and
  with dw_feas [0 0 0] (base 0 + raise loop) to pick the base for the direct form.
- 05:20 Phase B verification: default output bit-identical; the per-space rule raised w in
  non-owned directions (D drops, w rises); decision: registry D_d, w_d in non-owned directions,
  per-space numbers in owned directions only; invariant w_{k,d} <= w_d dropped (a lower own D
  may need a higher w). Sizes measured before the fix: self-sized decay slack -14% (S18),
  -23% (FO2s); like-sized +27-31% (the non-owned effect). Ladder starts after the rerun.
- 05:25 Gap study: G3 control reproduces (187,984 infeasible); harness 'auto' at w = 1 is the
  faces at offset 0 (210,340).
- 05:13 Litbuild finished (record copied to mixed_domain/litbuild_results.md): item 4 certified
  to 0.9875 d* under stock light (2.10M variables, 626-909 s per solve); its slot released.
- 05:12 G3 closed in both forms: direct form, faces at offset 1 (facesS 56,548; per-direction
  91,972) certify, offset 0 ('auto' at w = 1) 210,340, none/box/product2 fail; storage form:
  nothing at light under stock (faces0, box, corners, dw_feas 3, wcross 0 all infeasible), but
  the heavy LENGTHS as a row (Rheavy_u 188,320; Rheavy_s 189,216) certify at light size with no
  raise (the 2-D entry d3 row 1 = 2 instead of 1). Candidate tier change for Phase C: light 2-D
  d3 [2 2; 2 2]; and the 'auto' offset switch (offset 1 only when every w_k >= 2) costs 4x here.
- 05:35 WP4b (B, dw_feas 0, own-support sizing): S3 654,339 certified, S18(0.5) 400,899,
  S8 1,280,596, S26 1,167,107 certified; S17(0.5) UNKNOWN 1,309,239 (same size as like-sized).
  Cause: DEG omitted gives the 1-D blocks of a mixed PIE the 1-D entry [2 2 2] of lpi.var.deg,
  and the registry-wide slack rule carries that degree into every space. Executive edit (MMP):
  Q at st.lpi.var.deg{min(nv,end)} on every block (harness convention); WP4b reruns the mixed
  rows. Settings consequence: the direct form reads one entry of lpi.var.deg (the registry
  dimension); the per-block rule of lpivar_cdopvar stays for the synthesis LPIs.
- 05:47 G2 (S18 at 0.9 tau*): direct form, smallest closure = free operator at degree 2 (Q2)
  with facesS at dw 0: 235,488 (Q2 + 'auto' = faces0: 940,188; Q1 + dw 1 + facesS: 440k, with
  wcross 0: 374,827; faces at per-direction offset 2: 390,435; box 399,604; none out). Storage
  form: only Rlight2 (veryheavy lengths), 953,001 at dw_feas 0 (1.51M at dw_feas 1); Rheavy_s
  unknown at 1.73M; faces0/box/dw_feas 3 on Rlight0 infeasible. Pattern with G3: the faces at
  offset 1 (facesS) certify where 'auto' picks offset 0 at 4x the size; the degree of the free
  operator is a closing knob (S18 needs Q degree 2). Phase C candidates: (i) 'auto' = offset 1
  with the raise loop stepping to offset 0 before dw+1; (ii) lpi.var.deg 2-D entry 2 as a
  raise step or tier value (report already proposes 2 in heavy/veryheavy).
- 06:05 WP4b: with the explicit degree (registry entry) and base raise 0, the executive's
  direct form reproduces the harness bit for bit: S3 383,270; S18(0.5) 227,510; S17(0.5)
  435,350 certified; S8 492,048; S26 300,395 (all certified, 14-22 s). Base 1 (dw_feas 2-D
  entry) gives 1.18-1.19M on S3/S18 and refuses S17 at 5.2M. Executive edit (MMP): the 'Q'
  form starts at base 0 (own-sized slack; dw_feas/dw_P/dw_Q are allowances for like-sized
  slacks; dw_Q is read by H2_norm, Hinf_gain, Hinf_OFcontrol, so it was not changed).
- 06:05 Phase B verification 2 (owned directions only): the per-space degrees equal the
  registry rule on every case and target (S3, S17, S18, AGE, FO2s, heat2d); the degree
  variant is moot on this suite; only the per-space raise (dw_space) remains to measure.
- 06:00 G1 (S17 at 0.9 r2*): direct form closes with Q degree 2 + facesS at 290,588 (agent's
  interim note; table pending). Ladder on G4 (item 9).
- 06:20 WP4b final (scratchpad wp4b/WP4B_RESULTS.md): DEG omitted at base 0 is 1.7-3.9x the
  harness and leaves S17 UNKNOWN; base 1 costs 1.8-2.9x more and gains no verdict; explicit
  degree + base 0 = harness bit for bit (PDE33 190,900; S3 383,270; S18 227,510; S17 435,350;
  S8 492,048; S26 300,395; all certified at 0 raises). The repository executive now has both
  (explicit degree 05:35, base 0 06:12); a final rerun through the repository file as it
  stands belongs to the WP5 baseline rerun.
- 06:45-08:20 The three Opus agents were stopped by the 5-hour usage limit (reset 08:20); their
  last MATLAB scripts finished on their own (gap G4 rows, Phase B AGE/S18 rows). Resumed 08:25.
- 08:25 G1 (S17 at 0.9 r2*) closed in the direct form: Q degree 2 + facesS + wcross 0 at
  267,968 (facesS with the rule 290,588); at Q degree 1 no slack raise certifies (dw 1, 2: U).
  G4 from the logs: FO2r at 0.1 tau*: Q2 + facesS certified (351,677; wcross 0: 302,957), Q1
  infeasible at every slack raise; storage form: Rheavy_s certified at 847,204 (0 raises),
  Rlight2 at dw_feas 0 1,419,841, stock infeasible. FO2s at 0.95 tau*: Q2 + facesS certified
  (293,997; wcross 0: 248,497), Q1 unknown throughout; storage form: nothing (stock
  infeasible, Rlight2/Rheavy_s unknown). Pattern on G1, G2, G4: the free operator degree is
  the closing knob of the direct form, the faces at offset 1 the right term set; the slack
  weight raise closes nothing. Storage form: heavy lengths close G3 and FO2r.
- 08:30 Phase C launched (Opus): default-test candidates C0 (Q1, auto), C1 (Q1, facesS), C2
  (Q2, facesS), C3 (C2 + wcross 0), C4 (Q2, auto) on the suite, the literature set, the library
  2-D and 1-D controls, the 3-D heat; record phasec_results.md; recommendation for the tier
  default of the direct form and the first raise step.
- 08:55 Phase B (perspace_results, S18 at 0.9 tau*): the per-space rule changes no size or
  verdict (owned-direction degrees equal the registry's), but a PER-SPACE RAISE closes the
  direct-form gap: the weight of the 1-D space {s1} raised by 1 alone certifies at 325,275
  (faces at offset 0) and at 93,991 with facesS, against 440,023 for the uniform dw 1 with
  facesS and 215,188 for the Q-degree-2 closure; raising the 2-D space alone certifies at
  1,084,250. Storage form: no slack configuration up to 1.72M certifies (the storage is the
  knob, Rlight2 953,001). AGE: single space, nothing changes. Candidate executive ladder:
  raise the slack of the lowest-dimension spaces first, then uniform (needs dw_space in
  get_lift_degs, forwarded by lpi_ineq_sop; a WP5/WP6 item with a tier value of
  max_raise_feas 2-D = 2).
- 08:55 G5 (3-D heat 0.9 kappa* of the C-series, kappa* = 1.5 pi^2): not closed under the 3M cap:
  'P' faces0 dw_feas 2 cap_delta 0 infeasible at 1,306,897; the executive 'Q' at light 4,069,288
  refused (the 2.71M certificate of the record is the Cor. 35 bench LPI, not the executive).
  Asked for two more runs: 'Q exec' with facesS, and with Q degree 2 (size first). G6 running.
- 09:30 G5: 'Q exec' with the faces at offset 1 fits (658,536) but is INFEASIBLE (229 s): in
  three variables offset 0 is needed at w = 1 (as on 11/15); Q degree 2 run pending. G6 (D2):
  direct 6.8M even with offset-1 faces and wcross 0 (refused); 'P' raised solve (2.51M) killed
  at the watchdog; open (size).
- 09:40 Phase B done (mixed_domain/perspace_results.md): the per-space degree rule is inert
  (owned-direction degrees = registry); wcross 0 saves 12-28% on every multi-space case, no
  verdict change; dw_space +1 on the 1-D space closes S18 direct at 93,991; no per-space raise
  closes S17/FO2s (Q degree 2 does) or any 'P' gap (the storage does).
- 09:45 WP6 launched (Opus): port wcross, dw_space and the named set 'faces1' (faces at offset
  1 whatever w) into get_lift_degs; lpi_ineq_sop merges wcross/dw_space from the tier and accepts
  a per-dimension cell for psatz; no perspace port. Decision pending Phase C: light 2-D terms =
  'faces1' (4x smaller than 'auto' at w = 1 on every 2-D case so far) with 3-D kept at 'auto'
  (offset 0 needed); heavy 2-D free operator at degree 2 (closes S17, S18, FO2r, FO2s).
- 09:45 Phase C first rows: 1-D controls at Q degree 1 fail PDE2/3/5 in the direct test (Q
  degree 2 certifies all; the 1-D tier entry [2 2 2] is right); PDE33: C1 (Q1, faces1) 52,660
  vs C0 190,900 vs C2 210,361 vs C4 840,121, all certified.
- 09:55 Gap study complete (mixed_domain/gap_results.md); DECISIONS D66-D71 written.
- 09:57 WP6 verified through the repository files: defaults bit-identical; S18 at 0.9 tau*
  certified in the direct test at 93,991 (dw_space [1 0], 'faces1') and at 65,731 with
  wcross 0; 'faces1' = codes [0 3 4 5 6] offsets [0 1 1 1 1] in 2-D. Note for the tiers:
  dw_space needs exactly M entries (one per space), so it is an application override, not a
  tier value.
- 09:58 Phase C suite rows so far: C1 (Q1 + faces1) certifies every case C0 certifies (S1,
  S3, S4, S8) at a quarter of the size (S3 95,930 vs 383,270; S8 123,168 vs 492,048); Q
  degree 2 (C2/C3) lands just above the eta gate on S3/S4/S8. S6 and S7 at 0.9 are
  INFEASIBLE under every candidate (the direct form fails where the storage form certifies:
  a new gap class of the default test; not laddered yet).
- 09:58 WP5 launched (Opus): the six tier files on the report's layout (rows, wcross 0,
  psatz_offset [], lpi.var.deg 2-D entry 2 in heavy/veryheavy, settings_2d and the four
  executive properties out), lpisettings_sop, exec_sop_settings defaults, exec_tools_sop 2-D
  branch without settings_2d, guards in the other readers; verification script over the 1-D
  and 2-D rows, the 3-D heat and one gain. psatz stays 'auto' until Phase C's suite rows
  decide the 2-D entry ('faces1' candidate).
- 10:35 WP5 verification (interim, new tiers with wcross 0): 1-D rows identical to HEAD ('P')
  and to WP4b ('Q'); PDE33 'P' 104,305 (separable row; Rdef 103,969), 'Q' 190,900 identical;
  S3 'P' 275,005 certified, 'Q' 334,010 (wcross 0: -13%); S8 'P' 244,521, 'Q' 383,328; S18
  at 0.5 'P' infeasible 428,185 after one raise, 'Q' certified 201,050; S17 at 0.5 'P'
  unknown, 'Q' UNKNOWN after one raise at 676,410 (the wcross-0 exception of T5 with the
  faces at offset 0). Decision: lpi.ineq.psatz = {'auto','faces1','auto'} in every tier
  (offset-1 faces keep S17 and certify every 2-D case at about a quarter of the size);
  WP5 reruns the 2-D rows.
- 10:30 WP5 complete (D73): tiers on the report's layout, verified; 2-D rows 3-4x smaller with
  'faces1' (PDE33 'P' 26,545, S3 84,205/86,078, S8 'Q' 101,424, S26 'Q' 69,915, S17 at 0.5
  'Q' 98,410 certified); the one regression S18 at 0.5 'Q' needs one raise (374,827). Phase C
  suite: S2, S6, S7, S9 at 0.9 infeasible under every direct-form candidate (the storage form
  certifies them): the direct form is not uniformly stronger; S6/S7 ladder requested.
- 10:40 Phase C2 launched: optional second output [PDE, st] of the example files with the
  closing settings (cases 43, 57, 58, 66, new 72 = S17); WP7 launched: settings_2d_stock(level)
  as the explicit object for the demos that read settings_2d, help text.
- 11:05 WP7 done: settings_2d_stock.m (18 of 18 level/workspace comparisons identical to the
  HEAD tier files); heat2d_compare, set2d_deg, pielr_lpi read it; demo margins unchanged
  (eppos applies in every dimension); help text of four executives, lpisettings_stock_sop and
  the README updated. Report-edit agent launched (Sections 5, 6, 11, new gap-study section,
  summary; rebuild).
- 10:50 Phase C2 done (D75): [PDE, st] in cases 43, 57, 58, 66, 72 (new: S17), verified.
- 11:10 Report rebuilt clean (44 pages; proposal 26): Section 5 (wcross, dw_space, 'faces1',
  psatz cell, direct-form sizing and base, settings_2d_stock, the [PDE, st] mechanism and
  tab:app-settings), Section 6 (knobs, "Implemented 11/17/2026"), new Section 11 (the gap
  study, the Phase C sweep, the complementarity of the two forms), Section 12 open items.
  Phase C's literature rows: TL (two-layer heat, case 64) at 0.9 lam* certified by the
  offset-0 faces (1,152,400) and UNKNOWN under 'faces1' (288,400): with S18 at 0.5 tau* the
  second case the 2-D tier entry hands to the raise loop; LB infeasible under every candidate.

- 11:50 Phase C complete (D76, mixed_domain/phasec_results.md): default C1 (= the tiers after
  D72) certifies the same 17 of 29 cases as today's 'auto' at a median 0.27 of the size; the
  ladder C1 -> Q degree 2 -> offset-0 faces certifies 24 of 29 at half the summed cost of the
  present ladder; S2, S6, S7, S9, LB are infeasible in the direct form under every knob and
  certified by the storage form (S6/S7: nothing under the cap, incl. dw 1 and Q degree 2):
  a fallback to the storage form is the remaining proposal. 3-D: every candidate above the cap.

- 12:05 Report final: sopvar/settings_report_2026-11-16.pdf (45 pages) and the proposal (26)
  build clean; Section 11 "The gap study of 11/17/2026" written from the final records,
  Section 12 open items, Section 5.8 the [PDE, st] mechanism with its table; Section 6 the
  tiers as implemented. Run closed at 12:05 (the window ran to 16:30; every campaign and
  package finished). Agents: gap study, Phase B, Phase C, WP4b, WP5, WP6, WP7, Phase C2,
  the examples library, the literature build, two report passes.

## Outcome (12:05)

Gaps closed: G1 S17, G2 S18, G3 AGE, G4 FO2r/FO2s (details in gap_results.md and D70).
Open: G5 3-D heat and G6 D2 under the 3M cap; the direct form behind the storage form on
S2, S6, S7, S9, LB (and TL under 'faces1'); the raise ladder (free-operator degree, then term
offset 0, then slack weight, lowest-dimension spaces first) proposed, not implemented.

## Uncommitted changes (nothing committed; the maintainer decides)

Repository code (36 modified files, 1,765 insertions, 1,134 deletions, plus new files):
- executives_sopvar: PIETOOLS_stability_sop (direct 'Q' form: own-support slacks, Q at
  lpi.var.deg{min(nv,end)}, base raise 0), exec_tools_sop (row route in every dimension, 2-D
  branch without settings_2d), exec_sop_settings (executive properties), H2_norm, Hinf_gain,
  Hinf_estimator (help text).
- lpi_programming_sopvar: poscopvar (row reader), poscopvar_direct (psatz_rows_perdir, cap),
  lpi_ineq_sop (tier merges, psatz cell), misc/get_lift_degs (wcross, dw_space, 'faces1');
  staged moves poscopvar_lift, posmult_cdopvar from orphaned_file/private to misc.
- settings_sopvar: the six tier files on the report's layout, lpisettings_sop,
  lpisettings_stock_sop, README, new settings_2d_stock.m.
- PIETOOLS_examples: Examples_Library/Mixed_Domain_Examples (17 files, cases 56-72), case 43,
  the index (n_examples 72).
- PIETOOLS_demos: DEMO12 (row), DEMO3, DEMO5 (earlier today), cuadmm and lowrank demos
  (settings_2d_stock), cuadmm_outdir.
- Records: sopvar/psatz_faces_results (DECISIONS D66-D75, gap_results, perspace_results,
  litbuild_results, phasec_results, this log), the report and proposal sources and PDFs.

## Decisions for the maintainer

1. Commit scope and order (code, examples, records, report) and whether the git mv of
   poscopvar_lift/posmult_cdopvar stays.
2. The 2-D tier term entry 'faces1' (D72): a quarter of the size on every 2-D case and S17 kept,
   against S18 at 0.5 tau* (one raise) and TL at 0.9 lam* (UNKNOWN, raise not measured).
3. The raise ladder of the executives (free-operator degree, offset 0, slack weight; per-space
   first) and a combined default test (direct form, then storage form on failure): design
   decisions with measured costs in the report (Sections 11 and 12).
4. The 3-D cases: offset 0 needed; G5/G6 only above the 3M cap (Sol).
5. The example-file mechanism [PDE, st] and the new case 72 (S17) in the library.

## Sol jobs (16:50, after the run; maintainer's word)

- 1846967 g5_heat3q: 3-D heat at 0.9 kappa*, executive 'Q' at light (4.07M variables).
- 1846993 g6_d2q: D2 at 0.5 tau*, executive 'Q' at light with the 3-D term entry 'faces1'
  (6.8M variables).
Check: ssh sol squeue -u mpeet; logs /scratch/mpeet/pietools/gap1117/jobs/g5_heat3q_1846967.out
and g6_d2q_1846993.out ('@@ H3Q', '@@ D2Q' lines, '## ... DONE'); results in gap1117/results.
- 17:13 Sol job 1846967 (g5_heat3q) DONE: the 3-D heat at 0.9 kappa* CERTIFIED in the
  executive's direct form at light, 4,069,288 variables, m 41,164, 1,426 s, 39.7 GB (D78):
  G5 closed above the workstation cap. D2 job resubmitted as 1851319 (script fix).
- 17:35 origin/ndopvar fast-forwarded into the checkout (HEAD 33d30877; D79); pathdef
  refreshed for the parser reorganization; regression rerun of WP5 part a on the merged tree.
- 17:40 Regression on the merged tree (WP5 part a, refreshed pathdef): every row identical
  to the WP5 final state (1-D 'P' = HEAD, 'Q' = WP4b; PDE33 26,545/52,660; S3 84,205/86,078;
  S8 244,521/101,424; LRM 0.0459441): Declan's eight commits change no program or verdict.
- 17:31 Sol job 1851319 (D2, executive 'Q' light with 3-D 'faces1'): the light 3-D entry's
  cap_delta 0 capped the slack to 475,404 variables, INFEASIBLE in 126 s (22 GB); resubmitted
  uncapped (cap_delta NaN, the 6.8M program) as job 1851724 at 17:40.
- 18:45 Sol job 1851724 (D2 uncapped: 6,813,852 variables, m 163,847) OUT_OF_MEMORY at 200 GB
  after 19 min, inside sossolve, no MOSEK status (D80): the dense normal-equation matrix alone
  is 8 m^2 = 200 GiB; a rerun needs 320-480 GB and 4-41 h (600-7,500 CHE), the maintainer's
  decision. G6 stays open in both forms; the rungs between cap_delta 0 and NaN sized next.
- 18:55 D2 rungs sized on the workstation (sol1117/size/d2_size.log, size-only lpisolve
  shadow, no solve; light 'Q', 3-D terms 'faces1', wcross 0, Q degree 1): the slack cap
  cap_delta(3) = 0: 475,404 variables, m 33,025 (INFEASIBLE, Sol 1851319); 1: 2,721,692,
  m 82,551; 2: 5,441,948, m 132,572; NaN: 6,813,852, m 163,847 (OUT_OF_MEMORY, Sol 1851724;
  the local build reproduces the Sol sizes exactly). By the 8 m^2 + 30 GB model the cap-1
  rung needs about 85-95 GB and, from the 10/30 two-field solve (m 85,344, 300-315 s per
  iteration on 64 threads), 2-6 h on Sol (200-600 CHE at 128 GB); the cap-2 rung about
  170 GB and 8-25 h (1,000-3,600 CHE). Neither submitted.
- 19:02 The two 2-D closure knobs on D2, sized in the executive's direct form at light (3-D
  'faces1', wcross 0; sol1117/size/d2_run.log; none solved: every rung is above the
  workstation gates ndv 3M, m 36,000): slack cap 0 + dw_space [1 0] (the 2-D space raised)
  1,415,292 variables, m 42,222; cap 0 + dw_space [2 0] 4,028,172, m 54,363; cap 0 + Q degree
  2 798,476 variables but m 119,511 (the free-operator degree multiplies the rows); cap 0 +
  both 3,075,980, m 133,202; cap 1 + dw_space [1 0] 4,447,820, m 94,491; cap 1 + Q degree 2
  6,520,956, m 172,797. The one rung inside the variable cap (cap 0 + dw_space [1 0]) is a
  Sol job of the 3-D heat's size (m 41,164: 20 MOSEK iterations, 1,352 s, 39.7 GiB, about
  45 CHE), not a workstation run (64 GiB, 35 free, while the maintainer is online).
- 19:03 MOSEK statistics of the Sol logs: the 3-D heat (m 41,164; after presolve 36,233)
  20 iterations in 1,292 s, factor 2.78e8 nonzeros (2.2 GB), 6.34e12 flops per iteration;
  C5d (m 71,349; 71,187) 20 iterations in 17,470 s, factor 1.31e9 (10.5 GB), 1.09e14 flops;
  D2 cap 0 (m 33,025; 22,295) 14 iterations in 97 s, PRIMAL_INFEASIBLE_CER. The killed job
  1851724 printed no MOSEK banner: it died in the SeDuMi-to-MOSEK conversion or the task
  input, not in the factorization, so the 8 m^2 model of D80 is superseded: peak memory
  follows nnz(A) (calibration in sol1117/size/calib_nnz.log), the factor is a minor term.
- 19:06 Maintainer: "Submit the cap 0 + dw_space [1 0] rung to Sol". Job 1858048 (g6_d2dws):
  1,415,292 variables, m 42,222, 64 cores, 120 GB, 4 h (94 CHE/h, expected about 0.5 h);
  size guard at exactly 1,415,292 (D81). Check: ssh sol squeue -u mpeet; log
  gap1117/jobs/g6_d2dws_1858048.out ('@@ D2W', '## D2W DONE'); result results/d2dws_sol.mat.
- 19:17 D2 tau bisection (capped light 'Q', 3-D 'faces1'): infeasible at 0.25, 0.125, 0.0625
  tau*, unknown below (to 0.0078 tau*); certified at no delay (D82): the 3-D cap, not the
  margin, is the block. Calibration: nnz(A) heat3 2.67e7 (39.7 GiB), D2 cap0 2.60e6 (21.9 GiB),
  cap0+dws10 5.96e6, cap1 1.58e7, uncapped 4.42e7 (>200 GB): memory follows neither 8 m^2
  nor nnz alone; cause open.
- 19:22 Sol job 1858048 (cap 0 + dw_space [1 0]) INFEASIBLE with a certificate: 1,415,292
  variables, m 42,222, 498 s MOSEK, 30.7 GiB, about 14 CHE (D83). The per-space raise does
  not repair the capped 3-D slack; the next rung is cap_delta(3) = 1 (2.7M, m 82.5k; Sol).
- 19:29 Library pass under the six tiers started (D84): harness copy in
  testing_codes/private/baselines/2026-11-17_tiers/private, 162 catalog rows (cases 53-72
  added), one MATLAB process, part A light/heavy/veryheavy then part B extreme/stripped;
  custom = light (not run); logs scratchpad libpass/logs/pass_{A,B}.log; results_A/B.csv.
- 19:40 exec_tools_sop: an explicit settings_2d wins in 2-D (POS_ROW not wired when nv==2 and
  st has settings_2d), after the low-rank session measured the row overriding
  settings_2d_stock (rd2d-deg3-f010 Gram [25 424] vs HEAD [8 424]); tiers unchanged (D86).
- 20:35 Memory model corrected (D90, Sister session): MOSEK's dense Schur complement, 23-25
  bytes per m^2 (the two solved jobs), not the conversion (<= 10 GB for the killed program):
  the uncapped D2 needs ~600 GiB (highmem), the cap-1 rung ~170 GB and 7-8 h on 64 cores.
- 21:57 HEAD reference run of the library pass complete (worktree C:/wt1117 at 33d30877, the same
  harness copy, catalog, size gate and MOSEK): light 20:19-20:57 and heavy 20:57-21:57, 158 rows
  each (results_light/heavy.csv in the worktree harness folder). Cases 56-72 error on HEAD (their
  example files are not at 33d30877); 3-D case 53 light certifies on HEAD (907,729 variables,
  268 s), 54 and 55 refuse at the gate. Comparisons: scratchpad libpass/cmp_vs_head_lh.md (new
  light and heavy against HEAD), cmp_head_vs_1028.md (HEAD against the 10/28 reference).
- 22:00 Three-way attribution, light and heavy. 10/28 -> HEAD: every DDE and NDSDDF value
  unchanged except DDE01 heavy (1.9085 -> 1.8587, same size) and status flips at the eta gate;
  m drops on every row (the 11/10 equality-row fix). HEAD -> new tree: DDE light thresholds
  down on DDE03 (-2.5%), DDE04 (-1.8%), DDE05 (-0.3%), DDE07 (-1.7%), DDE08 (-5.3%), NDS05
  (-0.08%); heavy unchanged on these (DDE01 heavy up 8.6%, to the light value); 2-D: 33, 34,
  35, 37 heavy at 3.6x fewer variables with the same verdict (751,849 -> 210,361), 40 (both
  tiers), 43 (both), 38 heavy, 13.2 heavy, 10 light, 31 heavy, 48 heavy now certified; 39 light
  certified -> unknown (feasratio 7.14, 648,900 variables, m 30,721, against 933,696 / 18,411),
  36 light and 38 light infeasible -> unknown. Cases 53-56 light and heavy must be rerun (the
  rows predate the lb_build_system guard).
- 22:07 Fixed-delay A/B on both trees (scratchpad libpass/ab/ab_dde.m, logs ab_dde_{new,head}.log),
  light tier, at four delays HEAD certifies and the new tree refuses (DDE05 1.5707, DDE08
  1.3998, DDE03 1.69101, DDE07 1.35191). New tree: infeasible or unknown under every slack
  knob (the ladder's two raises at 3,895-15,580 variables, dD 2, dw 2, dw_feas 1); certified
  on all four with lpi.var.deg{1} = [4 2 2], the heavy entry, at 1,549-6,196 variables and no
  raise. HEAD: certified on all four after its slack raises (883-1,852 variables; the base
  program infeasible on three), its Q at the degrees of the storage (get_lpivar_degs_sop).
  Cause: the direct 'Q' form reads Q from lpi.var.deg (PIETOOLS_stability_sop, the 11/17/2026
  line), light's 1-D entry [2 2 2] is too small on these DDEs, and the raise ladder raises the
  slack, which the own-support rule already saturates, never Q. The 2-D analogue (light
  var.deg{2} = 1 against heavy 2) is being tested on PDE 39 (ab_2d.m).
- 22:35 2-D A/B on PDE 39 light (2D_NS_Antonelli; scratchpad libpass/ab/ab_2d.m, logs
  ab_2d_{new,head}.log), one solve per variant, no raise. HEAD: the base program certifies
  (933,696 variables, m 18,411, 306 s; residual 3.9e-6, eta 1.2e-10, OPT_STATUS 1.00). New
  tree, the tier (var.deg{2} = 1, 'faces1'): 648,900 variables, m 30,721, 266 s, MOSEK OPTIMAL,
  residual 4.3e-5, eta 2.4e-9, OPT_STATUS 7.14: 'optimal' by the executive's gate, 'unknown' by
  the harness's OPT_STATUS gate (|fr-1| <= 0.3); the pass row is this base program (no raise
  ran). 'auto' terms: 2,592,900 variables, 805 s, the same split (residual 4.1e-4, OPT_STATUS
  4.94): the term set is not the cause; faces1 is 4x cheaper. var.deg{2} = 2 (heavy): 3,736,449
  variables, m 83,708 (13.3M with 'auto'), above the gate: not solvable on the workstation.
  So the direct 'Q' form at light is a conditioning loss on this case (the storage form of HEAD
  is cleaner at 1.44x the variables), not a verdict loss. PDE 36 light likewise: new tree
  optimal after one raise (1,775,200 variables, residual 3.5e-6, eta 6.3e-9, OPT_STATUS 4.05)
  where HEAD is infeasible; PDE 38 light optimal at the base (245,200, OPT_STATUS -0.12) where
  HEAD is infeasible. Harness gate against the executive gate over the pass logs: 'unknown by
  the harness, optimal by the executive' on 13 light, 3 heavy, 1 veryheavy rows of the new tree
  (HEAD: 8 light, 8 heavy); the reverse never occurs.
- 22:38 Light with lpi.var.deg{1} = [4 2 2] on the ten 1-D bisection rows (scratchpad
  libpass/ab/ab_dde_fix.m, the harness bisection): identical to the new heavy tier in value and
  size on all ten (DDE01 2.0180859, DDE02 2.0839844, DDE03 1.7161724, DDE04 6.1665521, DDE05
  1.5707, DDE06 3, DDE07 1.3721821, DDE08 1.42229, NDS05 1.5784509, NDS08 0.027388574; last
  programs 1,712-24,784 variables), so in the 'Q' form of a 1-D PIE var.deg is the only tier
  knob that acts. Against HEAD light: equal or higher on eight (DDE01 +0.4%, DDE03 +0.1%, DDE07
  +0.06%, DDE08 +0.3%), lower on NDS05 (-0.08%) and NDS08 (-0.1%); bisection times 5-47 s against
  3-21 s. Ladder cost seen in the pass: veryheavy DDE17 (max_raise_feas 3 in 1-D) runs four
  solves per infeasible bisection step (591,950 -> 865,110 -> 1,192,770 -> 1,574,930 variables,
  about 450 s per step, every raise infeasible with a certificate, eta 2e-2 to 5e-2): the
  bisection takes about 50 min, three quarters of it raises. A raise does flip verdicts
  elsewhere (PDE 43 light: infeasible at dw 0, optimal at dw 1).
- 23:07 Maintainer: "pause any unstarted testing runs and resume in 1 hour": two pause markers in
  scratchpad matlab_slots (the launch protocol counts files), a background timer removes them at
  00:07. 23:29 "pause at the next boundary": PAUSE_sop in the harness folder (lb_run stops before
  the next row). 23:32 "pause the bisection loop now": the pass MATLAB (PID 70732, veryheavy
  NDS09, third bisection step) suspended in place with NtSuspendProcess (CPU time frozen at
  1213.5 min, 80 threads waiting, 8.2 GiB kept); sossolve sets no MOSEK wall-clock limit, so the
  solve continues on NtResumeProcess. The part A watchdog (lb_pass.sh, 21600 s from 19:29) kills
  the process at 01:29 regardless; after a 00:07 resume part A ends about 00:45.
- 00:36 (11/18) Maintainer: "Task your subagent with designing settings files which are example
  specific, work better than the default tiers for that application (e.g. lower cost for same
  accuracy), and are set inside the example file." One Opus agent launched: scope in order the
  1-D delay rows DDE03/04/05/07/08, NDS05/08 (and how the DDE/NDS libraries can return a
  per-case st), the 2-D examples 33-40 (PDE 39: a clean certificate by both gates), the
  mixed-domain cases with an st block (43, 57, 58, 66, 72) and the gate refusals (50, 61, 65,
  69, 70), then PDE 45-48; st blocks inside the example files (nargout > 1), own results files,
  second MATLAB slot only, sessions under 20 min; note
  sopvar/psatz_faces_results/example_settings/EXAMPLE_SETTINGS.md; about 3 h.
- 00:50 (11/18) Correction by the maintainer: "Q is a generalized storage where P form is a
  special case when Q = PT." The direct form 'Q' is not storage-free: its storage is T'Q, with
  the P form the special case Q = P T (storage T'P T). Where D92, D93 and ab_2d.m say "no
  storage" they mean that the storage ROW knob lpi.pos is not read by the 'Q' branch, which is
  a code fact, not a statement about the certificate. Consequence for the degree entry of Q:
  the Q family contains the P form of the same tier only when deg Q >= deg P + deg T per
  direction; at light's [2 2 2] with the storage of kernel degree 2 and a T of degree 1 on the
  DDE PIEs it does not (deg PT = 3), which is consistent with the measured losses and with
  [4 2 2] restoring them (not verified symbolically).
- 00:58 (11/18) Maintainer: "make the separation" (D95): sop.deg_Q, max_raise_Q, raise_Q for the
  direct stability form (PIETOOLS_stability_sop, exec_sop_settings, exec_tools_sop, the six
  tier files, README), tier values those of 11/17 (programs unchanged); verification ab_sep.m
  queued. 01:28 Maintainer: "Set light's deg_Q to [4 2 2] and rerun the DDE rows": light's
  deg_Q{1} = [4 2 2] (custom unchanged, no longer equal to light); the DDE and NDSDDF light rows
  rerun by the harness into results_ddeQ4.csv (results_A.csv is still being written by part A).
- 01:45 (11/18) Maintainer: "propose differentiations between Q and P form in all the tiers":
  sopvar/psatz_faces_results/q_p_differentiation_proposal.md (deg_Q by the rule deg Q >= deg P
  + deg T per tier, max_raise_Q, raise_Q 'deg' for the light tiers and 'both' for heavy and
  veryheavy, a form ladder 'Q' then 'P' in 2-D, a cap after two infeasibility certificates;
  four measurements before deciding). Also 01:30: part A killed by its watchdog at 21,600 s
  (the 43-minute suspension counted), during PDE 56 veryheavy; veryheavy rows 56-72 queued as
  a resume (lb_tiers_1117('veryheavy')) behind part B and the 53-56 rerun; 01:30 the cuADMM and
  low-rank sessions informed of the separation and of light's deg_Q.
- 01:50 (11/18) Maintainer: results and tests belong in the cleanup folder, not sopvar/: the
  untracked results, tests and reports of sopvar/ moved to cleaning/private/settings_1116/
  (psatz_faces_results/ with DECISIONS.md and this log, the two 11/15 md files, report/ with
  the 11/16 report and the proposal); the harness stays in testing_codes/private/baselines.
  Also: the veryheavy resume withdrawn ("postponed until the more urgent tasks are done");
  ab_sep verification passed (D95); the low-rank session reports every fingerprint identical
  (its builders pose 'P'/'Pd' and never the executive's 'Q' form).
- 02:00 (11/18) cuADMM session re-measured its suites on the deg_Q settings (nothing re-banked):
  only the three light 1-D 'Q' rows moved and now equal the heavy 1-D programs, all certified
  (stab_rd1 m 105 -> 123, nnz 1,284 -> 2,192; stab_tr1 m 100 -> 118; stab_wave1 m 418 -> 488,
  nnz 5,261 -> 8,978); every heavy, 'P'/'Pd'/'Qd', 2-D, objective and scaling row structurally
  identical to 11/17; its private 'Q' mirror (stab_mirror) now reads sop.deg_Q with the
  lpi.var.deg fallback; the bridge demos pass with the new fields; cuADMM's 1-D stability
  residual 1.9e-6 at 123 rows (2.9e-7 at 105 on 11/17; the (1+||b||) normalization gap).
- 02:10 (11/18) origin/ndopvar fast-forwarded into the checkout (D97): HEAD dd465ee5; DEMO12 and
  the examples library three-way merged with the uncommitted edits, verified identical; the
  running jobs (part B, the DDE rerun, the agent) use the examples library, whose content
  changed by one leading space only.
- 02:18 (11/18) results_ddeQ4.csv: the DDE and NDSDDF light rows at deg_Q [4 2 2] (D98): the five
  thresholds back at heavy's values, NDS11 certified, NDS08 -0.5% against light [2 2 2] (-0.1%
  against HEAD); 474 s for 66 rows.
- 02:20 Shivakumar's bounds as knobs (maintainer: "make the block-specific cap change. Then add the
  measurements to the queue"): poscopvar_direct degree field cap_inc (theta_k + s_k on the lower
  and upper basis operators only), get_lift_degs opts.Dmax (D = Dmax_d) and opts.cap_inc ('auto'
  = Dmax_d - 1), lpi_ineq_sop passes both from the tier (lpi.ineq) or the call; the six tiers
  carry Dmax false, cap_inc [] (programs unchanged); README bullet. Queue in the second slot:
  ab_slack check (the four reference rows under light as is and the three variants), then
  Dmax, capinc, both on DDE, NDSDDF and the 1-D PDE rows (logs scratchpad
  libpass/logs/ab_slack_<job>.log).
- 02:33 (11/18) Demo pass before the commit (verify_demos.sh, each demo its own MATLAB process):
  tiers ok; DEMO7 ok (35 s); DEMO3 FAILS at lpisolve(prob): sossolve refuses the empty solver
  '' every tier carries since the 11/17 layout (resolved by exec_sop_settings for the
  executives only) -> lpisolve now resolves '' as exec_sop_settings does (MMP, 11/18/2026);
  DEMO12 (Dan Braghini's 11/18 version, merged at 02:10) FAILS at its line 138, "for dw =
  0:settings.sop.max_raise": the tier's max_raise is a per-dimension vector ([1 1 0] since
  10/29/2026; the executives reduce it in exec_tools_sop), his inlined loop treats it as a
  scalar. Dan Braghini's code: reported, not patched (the same on origin dd465ee5; the author is Dan Braghini, DB, not Declan). The executive
  check of the pass used the interactive examples library (not batch-safe) and is rerun
  through the harness (ab_exec.m).
- 02:45 (11/18) Example-settings agent finished (D99): 14 st blocks, large savings at the library
  default parameters (PDE 33 at 409 variables against 52,660), DDE08 at 3.9144; in 1-D the
  form 'P' at light beats every form-'Q' setting by 1.4-4.6x in variables, in 2-D 'Q' wins.
  Correction: deg T reads 0 on the DDE rows (agent's routine), so the containment rule does
  not explain the light DDE loss; the minimal Q degrees measured ([3 2 2] to [5 2 2]) stand
  as facts without that explanation.
- 02:55 (11/18) Maintainer on a form field in the settings: no; the forms prove different forms
  of stability (PDE-level), the example file specifies the executive option, and the 'P' form
  is used whenever it certifies, being the stronger form (D100). The proposal note is amended.
- 03:00 (11/18) Maintainer: the auto-executor is NOT changed to pose 'P' first, and the 1-D rows are
  not rerun ("actually, lets not do that"); the open item of D100 is closed without action.
  Nothing was edited for it.
- 03:12 (11/18) COMMIT bb367fe1 (D101) after the demo pass: tiers ok, DEMO3 ok after the
  lpisolve fix, DEMO7 ok, executives on PDE 1/20/33 at light and heavy reproduce the pass;
  DEMO12 fails in Dan Braghini's loop (reported). Left for the later commit: the agent's 14 st blocks
  and its note, the 'P' executive option in the 1-D delay examples (D100), the live part-B
  results, the proposal note, these records.
- 03:10 ab_slack results (D102): the triangular cap -17% variables on the 1-D rows at light
  with 6 better / 8 worse rows; the Dmax variant as implemented lowers w (invalid), redone
  with w held.
- 03:20 (11/18) PUSHED on the maintainer's word ("push"): origin/ndopvar dd465ee5 -> bb367fe1
  (fast-forward). Later commit still pending: the agent's st blocks and note, the 'P' executive
  option in the 1-D delay examples, the corrected lift-ceiling results, part B's final
  results, the proposal note, these records.
- 03:25 (11/18) Attribution corrected: the DEMO12 rewrite (78ab4501, dd465ee5, markers "DB
  11/17/2026") is Dan Braghini's, not Declan's; the records above are amended (the pushed
  commit message bb367fe1 says "Declan's 11/18 version" and cannot be changed). Report for
  Dan: cleaning/private/settings_1116/DEMO12_report_1118.md.
- 03:35 (11/18) DEMO12 fixed on the maintainer's word, inline format kept: max_raise and dw_Q
  reduced to PIE.dim before Dan Braghini's raise loop (the executives' reduction of
  EXEC_TOOLS_SOP); the loop, its exit test and the slack options read the scalars. Full run
  (synthesis + PIESIM) in progress in a held slot (scratchpad verify/demo12.log).
- 03:30 (11/18) Part B ended (D104): extreme 158 rows, stripped 146 of 158 (watchdog at 7,200 s);
  the whole pass is in results_all_1118.csv (harness folder) with the light DDE/NDSDDF rows at
  deg_Q [4 2 2]; final_summary_1118.md the per-tier summary against 10/28. Gaps: veryheavy
  56-72, stripped's last 12 rows, cases 53-56 at light and heavy (all postponed).
- 03:40 (11/18) DEMO12, second defect: its R and V came from POSLPIVAR_SETTINGS_SOP (stock
  fields mult/int/joint), which the tiers no longer carry; replaced by POSCOPVAR with DEG
  omitted (the executives' POS_ROW). Rerun in progress.
- 03:42 (11/18) Callers of POSLPIVAR_SETTINGS_SOP left: in the toolbox only EXEC_TOOLS_SOP (the
  no-row fallback); in testing_codes/private the tests e2e_hinf_gain_1d, test_dpvar_hinf_chain_sop,
  test_getsol_sop, test_translators_sop and cx_PIE2PDEstability call it with whatever settings
  they build and fail with the tiers unless they pass stock-format fields (not run; off the path).
- 03:45 (11/18) DEMO12 runs (D105): synthesis 5.6253 certified at dw 0; analysis 5.6237 optimal
  with the witness at 5.6486 (omega 7.8e5) above both bounds, so the guard prints no bound;
  PIESIM completes. Two fixes, inline format kept; uncommitted.
- 03:50 (11/18) Corrected ceiling measured (D106): 2.42x variables for bounds within 0.1%, 6 better /
  4 worse by numerics; the cap -17% with 6 / 8; both 1.19x with 8 worse. Tiers keep both off.
  All background runs of the night are complete; the agent is done.
- 03:55 (11/18) DEMO12 fix committed (03a8a261) and pushed (D107).
- 04:00 (11/18) Example settings committed (eec4f07d) and pushed (D108).
- 04:05 (11/18) Maintainer: "Have your subagent complete its study": the example-settings agent
  resumed on top of eec4f07d for PDE 36, 45-48 and the size-gated 50, 61, 65, 69, 70 (both
  slots free; the form reported, not set; Dmax/cap_inc per example only; about 2 h).
- 04:20 (11/18) Studies handed to the Sister session (D110, HANDOFF_1118.md): A D2 and the 3-D
  cap, B tier degrees, C extreme, D the 'P' option mechanism, E the pass gaps, F the report pdf;
  its first deliverable is a plan. The Fable session keeps the running Opus agent and relays.
