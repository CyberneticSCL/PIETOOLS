# Audit of records made with substituted forms (11/21/2026)

Maintainer, 11/21/2026: "Yes, audit". Scope: records that report PIETOOLS_stability_sop in form 'Q' on a PIE in two
or more variables while the executive solved 'P' for it, and 'Qd', 'FQ' or gain 'Qd' on such PIEs before 905b0d26.
Method: git history of the executive; a read-only scan of every memory file and every markdown record under
cleaning/private/ and testing_codes/private/ by an Opus agent (about 100 files), with its main citations checked
here (D91 line, GAP_RUN_LOG 05:06 entry, EXAMPLE_SETTINGS.md lines 7 and 354).

## The windows

| Substitution | Where it applied | Window |
|---|---|---|
| 'Q' -> 'P', PIEs in two or more variables | main working tree | 10/28 04:12 (2e3fa94e; two variables only until 10/29, when the executives first took 3-D PIEs) until the direct 'Q' entered the working tree on 11/17 (WP4 programs by 05:06; the base-0 version, D68, at 06:12) |
| 'Q' -> 'P' | worktree C:/wt1117 at 33d30877, the HEAD reference of the 11/17 library pass | 11/17 20:19 – 22:35 (and any later run on that tree) |
| 'Qd' -> 'Pd', 'FQ' -> 'FP', gain 'Qd' -> 'Pd' | every tree | 10/28 until 905b0d26, 11/21 02:05 |

Main-tree records between 11/17 05:06 and bb367fe1 (11/18 03:11) are genuine direct-'Q' results: the 11/17 tier
pass, D99 and EXAMPLE_SETTINGS.md ran on "33d30877 plus the uncommitted 11/17 tier work".

## Result

Nearly every affected conclusion traces to one source: the HEAD reference run on C:/wt1117 (11/17 20:19–21:57)
and its ab_2d comparison (22:35), whose 2-D 'Q' column was 'P' but was read as the old storage-form 'Q'. Records
written inside the main-tree window state 'P', call the executive with an explicit 'P', build their own LPIs, or
ran through the qform shadow of the 11/16 studies (pf_setup.m), which lifted the substitution. No 'Qd', 'FQ' or
gain 'Qd' result on a PIE with two or more variables is mislabelled except the 11/20 'FQ' rows of library case 43
(PIETOOLS_PDE_Ex_Heat_Eq_w_Interior_Delay): the catalog has no dual row on such a PIE.

### Conclusions that change (class b)

| Record | When | Statement | Reading with 'P' on the HEAD side |
|---|---|---|---|
| settings_1116/psatz_faces_results/DECISIONS.md D91 (line 97); restated in mixed_domain/GAP_RUN_LOG.md 322–334, HANDOFF_1118.md 31–34, memory settings-report-1116.md, report libpass.tex | 11/17 22:00 | The working tree "certifies 40, 43 (both tiers), 38 heavy ... that HEAD did not, at 3.6x fewer variables on 33/34/35/37 heavy, and loses 39 light"; "53 light certifies on HEAD (907,729)" | A comparison of the direct 'Q' (working tree) with 'P' (HEAD), not of two versions of 'Q'. 40, 43 certify in direct 'Q' where 'P' did not; 39 light and 53 light certify in 'P'. Of the "2-D" list, 13.2, 10 and 31 are one-variable rows and 48 the estimator |
| DECISIONS.md D93 (line 99); copies in GAP_RUN_LOG 346–360, ab/README.md, libpass.tex | 11/17 22:35 | PDE 39 (PIETOOLS_PDE_Ex_2D_NS_Antonelli) light: "a conditioning loss of the direct 'Q' form, not a verdict loss"; "HEAD's storage-form program certifies cleanly (933,696 variables)" | HEAD's program was 'P'. The comparison is between forms, so "conditioning loss of the same test" does not follow; the direct 'Q' program (648,900 variables, MOSEK optimal, eta 2.4e-9, measure 7.14) is a candidate for the a posteriori check |
| psatz_faces_results/q_p_differentiation_proposal.md 17–65 | 11/18 01:45 | HEAD is "the old storage form"; "'P' was not measured"; proposes a 'Q' -> 'P' fallback ladder and to "measure the 2-D fallback on 36, 38, 39 light" | HEAD had measured 'P' light on all three. The ladder was withdrawn at 02:55 for other reasons (D100) |
| psatz_faces_results/example_settings/EXAMPLE_SETTINGS.md 354 | 11/18 02:00 | "The storage form of HEAD (933,696, certified) cannot be selected through st" – so no st for PDE 39 | That program was 'P', which the example can name (D100); 39 certified in 'P' on 11/20 (problem_cases_1120 §6) |
| example_settings_1119/PHASE_C.md 527–528; memory autonomous-run-1119.md 86 | 11/19 19:00 | PDE 39 certified "on HEAD 33d30877 by the storage-based 'Q'" | It was 'P' |

### Mislabelled, nothing drawn (class a)

- 2026-11-17_tiers/head_33d30877/cmp_vs_head_lh.md 122–163: the HEAD column ('P') beside the working tree (direct
  'Q') as one test – the data behind D91.
- head_33d30877/cmp_head_vs_1028.md 189–230, 260: both columns 'P'; line 260 labelled 'PIETOOLS_stability_sop/Q'
  (the label was set before the substitution).
- 2026-11-17_tiers/private/final_summary_1118.md 1–8: same/differ counts against 10/28, whose 2-D rows were 'P'.
- memory container-executives-sop.md 22–24 (10/28: 42 certified, 43 infeasible – 'P') and mixed-domain-heavy-gains.md
  36 (36, 38 "infeasible (1650 s)" – 'P').
- 2026-10-28_example_library sop_results.csv (25 two-variable 'Q' rows, 18 certified), 2026-11-02_opvar_phase0 F3
  passes (8 such rows each), head_33d30877 results (11 light, 6 heavy): every two-variable 'Q' row is a 'P' result.
- The 11/20 'FQ' rows of case 43: corrected in problem_cases_1120; memory autonomous-run-1120.md 55–56 corrected
  here.

### Checked and unaffected (class c)

About 65 record files and 31 memory files: 'P' stated or requested explicitly (D49, the 3-D heat of 11/16, ran
'stabP' = form 'P'), custom LPIs (heatNd_lpi, the pf_* direct test), the qform shadow, one-variable PIEs (41, 44,
44min), or after 905b0d26.

## Effect on the settings

No example setting rests on these records: the 11/19 search took its two- and three-variable choices from the
11/18 studies on the main tree. The decisions above bear on three cases, all in the rerun (settings_rerun_1121,
starts_P1.csv): 39 (an alternative start in 'P' with the stock 2-D storage, besides 'Q'), 53 (alternatives in
'P' at light and heavy) and 36.

## Corrections made

A dated correction line was appended to DECISIONS.md, GAP_RUN_LOG.md, q_p_differentiation_proposal.md,
EXAMPLE_SETTINGS.md and PHASE_C.md (history left as written), and the memory files settings-report-1116,
autonomous-run-1119, autonomous-run-1120, container-executives-sop, mixed-domain-heavy-gains and
stability-form-swap-fix were corrected or annotated.

Not read: the settings_1116/report .tex and .pdf files beyond libpass.tex lines 23 and 51; scratchpad copies of
the libpass, wp4 and gapstudy notes; cleaning/private/PIETOOLS_clone_wcomments (a code clone).
