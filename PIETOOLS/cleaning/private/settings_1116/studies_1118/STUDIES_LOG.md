# Settings studies of 11/18/2026: run log (Sister session)

Plan: ../PLAN_1118_studies.md, approved by the maintainer on 11/18/2026. Decisions given with the
approval:
- the settings hook in a copy of the harness;
- a change of verdict counts only where both gates agree;
- D-a for the form option;
- the report written by this session;
- tightness data before any per-dimension dD.

Maintainer, 11/18 04:31: "run autonomously for the next 10 hours, following the proposed plan and
obtaining as many numerical testing results on example files as possible".

Harness: testing_codes/private/baselines/2026-11-18_studies/private, a copy of the 11/17 library
harness. Its additions:
- st_run.m (one configuration per tag) and the hooks in lb_exec.m (LB_ST_HOOK, LB_FORM_MAP);
- LB_ROW_KEY and opts.catalog in lb_run.m;
- the executive status, eta and passes in the note of every row;
- cmp_runs.pl, the comparison under the both-gates rule.

Results: results_<tag>.csv in the harness folder; MATLAB logs: logs/<tag>[_<attempt>].log here.

Rules kept:
- one MATLAB slot through the lock files of the Fable scratchpad matlab_slots (sister_<tag>.lock);
  a second queue only while no other session holds a lock;
- MATLAB sessions near 20 min (stop time 15 min; lb_run resumes);
- no Sol job, no commit, no push without the word of the maintainer.

## Schedule and estimates (from the unit costs of the 11/17 pass)

| block | content | estimate |
|---|---|---|
| Q1 | G46P; L1ctl; B2dD0; H1ctl; B3q522; B3q622; D1P; D2st; E1str; E2l; E2h | 1.5–2 h |
| Q2 | B1 tightness T0–T3 on 7 systems (28 bisections) | 2–4 h |
| Q3 | B1 candidates on the 22 2-D library rows and the 16 mixed rows | about 2.3 h |
| Q4 | C extreme, 1-D candidates then the best on 2-D | 0.5 h |
| Q5 | B4 veryheavy deg_Q{1} [8 6 6] on the ten 1-D bisection rows; E veryheavy 56–72 | 1–2 h |
| A2 | D2 size sweep, after the d2_size2 numbers of the Fable session | 15 min |
| D-a | stability_form patch after D1P/D2st, verified through PIETOOLS_DDE | 30 min |

## Log
- **S1 (04:37, G46P).** PDE 46, the 2-D L2 gain with exact value 0.17110, at light in form 'P' of
  PIETOOLS_Hinf_gain_sop.
  - Result: certified gamma = 0.17173673 (+0.37%) at 200,986 variables, m 3,983, 24.5 s; MOSEK
    OPTIMAL, optimality measure 1.00, eta 3.1e-10, one pass.
  - The same row in form 'Q' (the form of the auto-executor and the pass) is UNKNOWN in every tier:
    light gamma 812.6 at 351,307 variables (numerr 2), heavy 8.51 at 1,399,419, the same on HEAD and
    on 10/28.
  - So the 2-D gain failure of the library is a failure of the 'Q' form of the gain executive, not
    of the tiers. PDE 45 and 47 in form 'P' are queued (G45P, G47P).
- **S2 (04:49, L1ctl).** Light on the 124 one-variable rows (PDE, DDE, NDSDDF), 11 min. This is the
  control of B2 and D.
  - Check of the harness copy: the 66 DDE and NDSDDF rows reproduce results_ddeQ4.csv of 11/18 02:18
    bit for bit in result, variables and solves (897,139 variables, 221 solves).
  - The only differences are the gate splits already known (D93): DDE12, the dual of DDE26 H2, and
    the dual of NDS11 are unknown by the harness and optimal by the executive (eta 7e-10 to 8e-9).
- **S3 (same run).** The 58 one-variable PDE rows at light with deg_Q{1} [4 2 2] against the pass
  (deg_Q [2 2 2]):
  - every verdict and every bisection threshold equal;
  - PDE 68 (Batty 2016) certified by both gates (12,564 variables, eta 4.2e-9) where [2 2 2] was
    unknown;
  - summed variables 443,760 -> 565,254 (1.27x).
  - Four rows split the gates as before (the duals of PDE 20 and 49; PDE 22, 27).
- **S4 (04:58, reported by the Fable session).** Second pass of its example-settings agent;
  EXAMPLE_SETTINGS.md Section 11 gives both gates per row.
  - PDE 45, 46, 47 at light in form 'P' certify 0.111982, 0.171737, 0.208535 against the exact
    0.11101, 0.17110, 0.20668 (+0.88%, +0.37%, +0.90%), at 200,986 variables each (fr of PDE 47
    0.891). In form 'Q' none certifies (numerr 2; gamma between 672 and 90,621).
  - PDE 48 certifies at heavy only (768,159 variables).
  - PDE 36, 50, 65 certify nothing. On PDE 65, 'Q' with dD 0 fits the size gate at 1,585,644
    variables, but MOSEK needs about 110 s per iteration and was stopped at iteration 5.
  - PDE 61 at light with dD 0 is certified at 197,415 variables, where light exceeds the gate
    (2,430,327).
  - The triangular cap certified nothing on these rows.
  - The example-settings agent has finished; from 04:58 it holds no MATLAB slot.
- **04:58.** The maintainer has been offline since 04:55 (Fable session). The 20-minute session rule
  applies while he is online only, so st_run reads its session length from ST_STOPMIN in the harness
  folder: 60 min from now, 15 again when he is back. Commits, pushes and Sol jobs above about 10 CHE
  are held until he returns and recorded here. Diagnoses go to the Fable session (its offer of
  04:58).
- **S5 (05:03, B2dD0).** Light with lpi.ineq.dD = 0 on the 124 one-variable rows, against the light
  control S2 (studies_1118/B2_dD0.md).
  - Size: summed variables 1,462,393 -> 1,037,045 (0.71x), seconds 654 -> 527, solves equal.
  - Lost by both gates (6): PDE 7 dual (certified -> infeasible); PDE 29, 31, 51, 52 (certified ->
    unknown, bounds 4,012 to 13,014); PDE 30 (certified -> inaccurate, 0.6499 -> 1.1253).
  - Gained by both gates: none. Gates split differently on 9 rows (listed in the file).
  - Loser bounds where both certify: DDE13 1.3653 -> 1.4470 (+6.0%), DDE18 0.3779 -> 0.4152 (+9.9%),
    PDE 20 +0.41%, DDE11 +0.08%; the other certified gains within 1e-4 relative.
  - Bisection thresholds: DDE01 +0.10%, NDS08 +0.39%, PDE 5 -0.10%, the rest equal.
  - Verdict under the adoption rule: dD 0 is not a one-variable tier change. It loses six
    certificates, and the smaller lift loosens two gains by 6 to 10%.
- **S6 (05:15-05:52, H1ctl, B3q522, B3q622).** Heavy on the 124 one-variable rows (B3_q522.md,
  B3_q622.md).
  - Control: H1ctl reproduces the heavy rows of the 11/17 pass bit for bit (1,918,160 variables, 332
    solves). Against light it has equal bisection thresholds at 1.31x the variables (D94).
  - deg_Q{1} [5 2 2]: no certificate gained or lost by both gates; NDS11 dual becomes unknown by the
    harness while the executive stays optimal. Thresholds DDE01 2.0181 -> 2.0360 (+0.89%), DDE03
    1.7162 -> 1.71785 (the documented value), NDS05 +0.08%, NDS08 +0.49%. Cost: variables 1.36x,
    seconds 1.43x.
  - deg_Q{1} [6 2 2]: the same certificates, and the same four thresholds raised (DDE01 +0.69%).
    DDE08 rises 1.42229 -> 3.179895 (+124%; documented 3.9-3.95, veryheavy [6 6 6] 3.8957). Cost:
    variables 1.50x, seconds 1.49x.
  - Under the adoption rule both keep every certificate and certify larger parameter ranges.
    [6 2 2] is the heavy candidate for the maintainer: DDE03 at its documented value and DDE08 at
    2.2 times the heavy threshold, for 1.5 times the variables of the one-variable rows. A tier change
    waits for his word.
- **05:55. Study A changes course.** The second d2_size2 run of the Fable session ran without its
  size-only flag and solved the cap-1 rung of D2 at 0.5 tau* (2,721,692 variables, m 82,551) for 4
  MOSEK iterations before its watchdog.
  - The factorization is sparse (dense dim 3, 4.39e13 flops), about 190 s per iteration within 64 GB,
    PFEAS 4.0 -> 1.8e-2 (d2/d2_size2_cap1_mosek_partial.log).
  - So the 23-25 B/m^2 model of D90 does not govern this rung, and a full solve fits the workstation
    in an estimated 2-3 h.
  - Scheduled as A3cap1 (st_d2solve.m), an exclusive run once queue 1 ends.
  - Runners from here (st_qrun2.sh): A (qA: gains, A2 sizes, tightness, extreme, veryheavy delay
    rows), B2 (second queue: the B1 library candidates and the veryheavy gap rows, kept out while an
    exclusive run waits or runs), C (exclusive: A3).
- **06:05. A3 target changed.** D2 sizes from the Fable session (d2/NOTES.md, d2_size2.log), size
  only, at 0.5 tau*:
  - cap 1: 2,721,692 variables, m 82,551; with the triangular cap 1,394,588, m 63,066;
  - cap 0, with or without the triangular cap: 475,404, m 33,025;
  - cap 2 with the triangular cap: 2,027,420, m 84,890;
  - uncapped with the triangular cap: 2,134,428, m 90,777 (uncapped alone 6,813,852, m 163,847).
  - With the triangular cap the negativity slack is 735,840 variables at every cap.

  By Lemma 1.2 of the Newton polytope note the triangular cap keeps every certificate on a purely
  integral target, and D2 has no ODE. So 'uncapped + triangular cap' decides what the uncapped light
  direct form certifies:
  - if it fails, no cap rung certifies D2 at the light degrees;
  - if it certifies, cap 1 + triangular cap (m 63,066) is the next and cheaper rung.

  The exclusive run is now A3unc. A memory watchdog (st_memwatch.sh) stops it above 56 GB. Estimate:
  m 90,777 at about 190-250 s per iteration (sparse factor, as on the cap-1 control) for 30-50
  iterations, 2-3.5 h.
- **S7 (05:53-05:58, D1P, D2st).** Study D, the stability form of the delay examples (D1_formP.md).
  - **D1P:** light in form 'P' ('Pd' for the duals) on the 66 DDE and NDSDDF rows, against light 'Q'
    (L1ctl).
    - No certificate lost by both gates. The dual of NDS11 is certified by both gates in 'P'; in 'Q'
      the harness gate refused it while the executive certified it.
    - Thresholds: DDE01 2.0181 -> 2.04 (+1.09%), DDE03 1.7162 -> 1.71785 (documented), NDS05
      +0.08%, NDS08 +0.49%; the rest equal.
    - Cost: variables 897,139 -> 509,271 (0.57x), seconds 503 -> 218.
  - **D2st:** the st blocks committed in eec4f07d, in 'Q'.

    | row | block st in 'Q' | light 'P' |
    |---|---|---|
    | DDE03 | 1.71785 at 3,337 variables | 1.71785 at 1,220 |
    | DDE04 | 6.1665521 at 6,184 | 6.1665521 at 2,604 |
    | DDE05 | 1.5707 at 1,549 | 1.5707 at 1,107 |
    | DDE07 | 1.3721821 at 6,196 | 1.3721821 at 1,356 |
    | DDE08 ([8 6 6]) | 3.9144165 at 30,356 | 1.42229 at 1,464 |

  - Per-row decision: 'P' for every stability test of the DDE and NDSDDF libraries, except Example E
    (DDE08), which keeps 'Q' with its st.
- **D-a implemented (06:05, uncommitted).**
  - PIETOOLS_DDE.m and PIETOOLS_DDF.m set stability_form = 'P' in their settings section and call
    PIETOOLS_stability_sop with stability_form ([stability_form 'd'] for the dual).
  - The Example E block of examples_DDE_library_PIETOOLS.m sets stability_form = 'Q' (commented, like
    its st).
  - PIETOOLS_auto_execute_sop reads stability_form, with 'Q' when absent, so its default path is
    unchanged.
  - Header entries, column-77 markers and (was) lines in all four files.
  - A deviation from the plan text: the plan assumed that the DDE entry script runs the
    auto-executor. It calls the executives itself, so the variable is read in the two entry scripts
    as well. These scripts set 'P' as the form of their examples, and the Example E block overrides
    it.
  - Verification queued as tag VDA: two copies of PIETOOLS_DDE.m with Example B and Example E
    uncommented, simulation off.
- **S8 (05:58-06:08, E1str).** Library gap: stripped on PDE 61-72 (results_E1str.csv, seeded with
  the extreme rows for the 2-D projection).
  - Certified: 62 (86,078 variables), 63 (165,182), 64 (288,400), 66 (56,548), 67 (18,650).
  - Gate split: 68 (6,714; harness fr 1.45, executive optimal).
  - Unknown: 71 (the 2-D gain in 'Q', numerr 2, gamma 1,138) and 72 (S17, numerr 2).
  - Refused by the size gate: 61 (2,430,327), 65 (2,845,420; m 67,500), 69 (4,453,560), 70
    (2,014,033).
  - The same pattern as light on these rows.
- **S9 (06:25-06:28, A2size).** D2 at 0.5 tau*, light 'Q', 3-D 'faces1', size only (A2size_1.log).

  | configuration | variables | m |
  |---|---|---|
  | cap 0 | 475,404 | 33,025 |
  | cap 0, dD 0 | 155,004 | 23,286 |
  | cap 1 | 2,721,692 | 82,551 |
  | cap 1, triangular cap | 1,394,588 | 63,066 |
  | cap 1, dD 0 | 655,708 | 35,881 |
  | cap 1, deg_Q{3} 0 | 764,582 | 42,665 |
  | cap 1, triangular cap, dD 0 | 387,420 | 29,348 |
  | cap 1 or 2, dD 0, deg_Q{3} 0 | 83,574 | 8,733 |
  | uncapped, triangular cap, dD 0, deg_Q{3} 0 | 102,582 | 10,569 |
  | cap 1, 3-D 'auto' (offset-0 faces) | 16,274,652 | 144,751 |

  The lift raise and the degree of Q in three variables decide the size of D2. Solves of every rung
  under m 45,000 are queued (A3s1-A3s6) as cheap direct tests at 0.5 tau*.
- **06:30. A3unc started** (uncapped + triangular cap, 2,134,428 variables, m 90,777, exclusive), on
  its own runner (qC2) to avoid waiting for queue 1. The B1c1 session of runner B2 was stopped to free
  the slot (about 6 rows done; it resumes from its results file after A3).
- **S10 (05:49-06:25, B1c0).** Light on the 30 catalog rows in two variables (2-D library and mixed
  rows), the control of B1. It reproduces the light rows of the 11/17 pass bit for bit (7,338,777
  variables, 24 solves).
  - Gate splits (harness unknown, executive optimal): PDE 36 (eta 6.3e-9), 38 (1.9e-9), 39 (2.4e-9),
    48 (3.4e-10).
  - Certified: 33, 34, 35, 37, 40, 42, 43, 57, 59, 60, 62, 63, 64, 66.
  - Unknown: 45, 46, 47 (gains in 'Q'), 71, 72.
  - Infeasible: 58.
  - Size gate: 50 (and its dual), 61, 65, 69, 70.
- **S11 (06:30-06:34, A3unc).** D2 at 0.5 tau*, uncapped with the triangular cap (2,134,428
  variables, m 90,777), stopped by the memory watchdog at 57.2 GB, 4 min 20 s after start.
  - Private bytes of the session: 24.0 GB at the build; 17.4 GB at the conversion; then 30.7, 37.3,
    43.9 and 57.2 GB at 30 s intervals during the MOSEK setup, before the first iteration.
  - The sparse factorization of the cap-1 control does not keep this rung under the workstation
    limit, so it is a Sol job (held: above about 10 CHE, needs the maintainer).
  - The cheap rungs of S9 (m 8,733-42,665) run in its place: A3s3 started 06:34.
- **06:40. A3cap1ci queued (exclusive, watchdog 56 GB).** D2 at 0.5 tau*, cap 1 with the triangular
  cap: 1,394,588 variables, m 63,066. By Lemma 1.2 it carries the certificates of the cap-1 program.
  - Memory (Fable session, from the MOSEK statistics of the cap-1 control): factor 1.41e9 nonzeros,
    about 11.3 GB, peak below about 45 GB; a fully dense factor at m 63,066 is 15.9 GB.
  - Estimate: setup about 5 min, then about 100-190 s per iteration for 30-50 iterations, so 1-2.5 h.
  - The dD-0 rungs are a different family: a certificate there is sound, but an infeasibility does not
    speak for the rule's program (D113).
- **S12 (06:34-06:36, A3s3).** D2 at 0.5 tau*, cap 1, triangular cap, dD 0 (387,420 variables, m
  29,348): PRIMAL_INFEASIBLE with a certificate after 111 s of MOSEK (126 s in all). dD 0 lowers the
  lift below the rule, so the infeasibility does not decide the rule's program; it closes this
  cheaper family at 0.5 tau*.
- **06:37. A3cap1ci started** (exclusive, watchdog 56 GB).
- **S13 (06:08-06:37, E2l).** Library gap: light on PDE 53-56 (three variables). Every row is refused
  by the size gate:

  | row | variables | m | build |
  |---|---|---|---|
  | 53, the 3-D heat | 4,069,288 | 41,298 | 67 s |
  | 54 | 24,921,000 | 120,749 | 273 s |
  | 55 | 99,684,000 | 493,339 | 1,306 s |
  | 56, D2 with the 3-D 'auto' terms of the tier | 2,524,300 | 51,558 | 44 s |

  PDE 53 is the program Sol certified in D78. Heavy on the same rows (E2h) was started at 06:37
  beside the exclusive D2 solve and stopped after 28 s: its rows are larger than light, so they would
  only rebuild programs for the gate to refuse, at up to an hour each. Not rerun.
- **S14 (06:37, VDA).** D-a verified through copies of the patched PIETOOLS_DDE.m (vda/vda_B.m,
  vda_E.m; one library block uncommented each, simulation off, executive verbose).
  - Example B of [1] (DDE03) at tau 1.5: stability_form 'P' from PIETOOLS_DDE. 'P' is optimal at 404
    variables, 44 rows, eta 7.0e-12; 'Pd' optimal at 404, eta 1.2e-11; spectrum max Re -0.1182.
  - Example E of [2] (DDE08) at tau 1.1: stability_form 'Q' from its block, with the block st [8 6 6].
    'Q' is optimal at 30,356 variables, 756 rows, eta 3.3e-11; 'Qd' optimal at 733, eta 1.8e-11;
    spectrum max Re -0.0403.
- **S15 (06:38-06:39, A3s2, A3s1, A3s6).** D2 at 0.5 tau*, the smallest rungs of S9; all three are
  PRIMAL_INFEASIBLE with a certificate:
  - cap 1, dD 0, deg_Q{3} 0: 83,574 variables, m 8,733, in 9 s;
  - uncapped, triangular cap, dD 0, deg_Q{3} 0: 102,582, m 10,569, in 11 s;
  - cap 0, dD 0: 155,004, m 23,286.

  With S12 (cap 1, triangular cap, dD 0: infeasible), every dD-0 rung tried fails at 0.5 tau*.
- **S16 (06:39-06:44, GallQ, GallP).** The 21 H-infinity gain rows of the catalog at light, form 'Q'
  (the form of the auto-executor) against form 'P' ('Pd' dual) of PIETOOLS_Hinf_gain_sop
  (G_gainform.md, G_gainform_table.md). Sent to the Fable session for the diagnosis of 'Q'.
  - 'P' certifies all 21 by both gates; 'Q' certifies 14.
  - Gained by both gates: PDE 45, 46, 47 (0.111982, 0.171737, 0.208535 against the exact 0.11101,
    0.17110, 0.20668) and PDE 71 (0.045944081, the LRM record; 'Q' 169.9).
  - Gate splits in 'Q' that 'P' certifies: DDE12, PDE 20 dual, PDE 22.
  - Bounds: 'P' equal or tighter on every row certified by both forms except PDE 24 (+0.01%). PDE 20
    0.50081 -> 0.50000019 (exact 0.5), DDE11 -0.37%, PDE 23 -0.06%, PDE 26 -0.07%.
  - Cost: variables 1,451,265 -> 860,211 (0.59x), seconds 200 -> 79.
  - Whether the auto-executor and the entry scripts should pass 'P' for gains is the maintainer's
    decision.
- **S17 (06:37-07:06, A3cap1ci).** D2 at 0.5 tau*, light direct form, cap 1 with the triangular cap
  (1,394,588 variables, m 63,066): PRIMAL_INFEASIBLE with a certificate.
  - MOSEK: 1,711 s, about 67 s per iteration on 24 threads, 54,381 constraints after presolve, factor
    7.14e8 nonzeros (dense dim 0, 1.45e13 flops); peak 38.0 GB.
  - By Lemma 1.2 the triangular cap keeps the certificates of the cap-1 program, so the cap-1 rung
    fails at 0.5 tau*.
  - With D82 (cap 0: certified only at no delay) and D83 (cap 0 + dw_space [1 0]: infeasible), every
    capped light direct form tried is infeasible at 0.5 tau*.
  - Open: the uncapped program, equal in certificates to uncapped + triangular cap (2,134,428
    variables, m 90,777).
  - Its estimate from this run is about 2.1x the memory (scaling m^2: about 80 GB) and 3x the time
    (m^3: about 1.5 h on 24 threads). That is a Sol job of about 100-200 CHE, held for the
    maintainer.
- **S18 (07:06-07:11, A3s4).** D2 at 0.5 tau*, cap 1 with dD 0 (655,708 variables, m 35,881):
  PRIMAL_INFEASIBLE with a certificate, 296 s.
- **S19 (06:44-07:27, B1 tightness, first three systems).** Harness bisection toward the known
  margin (bracket [0.5, 1.1] x margin, tolerance 1e-3 of the margin; tight_table.sh). Variables are
  those of the last certified program, after raises.

  | system | T0 light | T1 dD 0 | T2 dD 0, Q from 0 by degree | T3 dD 1, Q from 0 by degree |
  |---|---|---|---|---|
  | KISS (34), lam < 2 pi^2 | 0.9998, 522,100 var, 353 s | 0.9951, 120,740, 64 s | 0.9998, 522,441, 330 s | 0.9998, 210,361, 253 s |
  | parabolic (37), a < pi^2/2 | 0.9998, 72,100, 246 s | 0.9998, 120,740, 47 s | 0.9998, 522,441, 283 s | 0.9998, 72,100, 137 s |
  | AGE (66), b < 2.302685 | 0.9881, 401,572, 251 s | 0.9869, 120,740, 56 s | 0.9951, 187,641, 98 s | 0.9951, 210,361, 193 s |

  Reading so far:
  - dD 0 loses reach on KISS (-0.47%) and AGE (-0.12%) at 4-5x lower time.
  - The degree raise of Q (T2, T3) reaches further than light on AGE (0.9951 against 0.9881) and
    ties elsewhere.
  - T3 matches or beats light everywhere at 0.4-1.0x its variables at the reach and 0.56-0.77x its
    time.
- **S20 (07:36, A3P1).** D2 at 0.5 tau*, storage form 'P' at light, uncapped with the triangular cap
  (9,108 variables, m 12,826; sized by the Fable session, d2/d2_size3.log): PRIMAL_INFEASIBLE with a
  certificate after a 1.7 s solve. The cap is taken from the like-operator R, not from the target,
  so Lemma 1.2 gives no guarantee and the infeasibility is inconclusive.
- **Study A closed for the run (07:37).** Every workstation-feasible rung at 0.5 tau* is infeasible:
  - cap 0 (D82: certified only at no delay);
  - cap 0 with dw_space [1 0] (D83);
  - cap 1, through cap 1 + triangular cap (D119);
  - every dD-0 rung (S12, S15, S18; not conclusive for the rule);
  - 'P' light with the triangular cap (not conclusive).

  The decisive programs are the uncapped ones, both above the workstation:
  - 'Q' uncapped + triangular cap: 2,134,428 variables, m 90,777; passed 57 GB in setup; estimated
    80 GB, 1.5 h;
  - 'P' light uncapped: 3,893,508 variables, m 95,216.

  Both are Sol jobs of about 100-200 CHE, held for the maintainer.
- **S21 (07:27-08:00, B1 tightness, S18 and S17; also A3s5).**

  | system | T0 light | T1 dD 0 | T2 dD 0, Q from 0 by degree | T3 dD 1, Q from 0 by degree |
  |---|---|---|---|---|
  | S18 (43), tau < 1.034726 | 0.9893, 374,827 var, 343 s | fails at 0.5 (infeasible, 115,647) | 0.9400, 196,513, 157 s | 0.9916, 215,188, 369 s |
  | S17 (72), r2 < pi^2 | 0.6781, 98,410, 433 s | fails at 0.5 (numerr 2, 173,990) | fails at 0.5 (harness fr 0.51; executive optimal) | 0.7268, 267,968, 414 s |

  - On S17 the harness gate bounds the reach: most failed steps are UNKNOWN, not infeasible. The gap
    study certified 0.9 r2* at the same 267,968 variables by the executive gate (D70).
  - dD 0 fails S18 and S17 at half the margin. Per-dimension dD = 0 in two variables therefore costs
    certificates; the maintainer's tightness concern holds.
  - T3 (light lift, Q of the direct form starting at degree 0 and raised by degree up to 2) reaches
    further than light on S18 and S17.
  - A3s5 (D2 at 0.5 tau*, cap 1, deg_Q{3} 0; 764,582 variables, m 42,665): infeasible, 881 s.
- **S22 (08:04-08:18, B1 tightness, FO2s, PDE 57, tau < 2.881401).**

  | candidate | reach / margin | variables | time |
  |---|---|---|---|
  | T0 light | 0.7420 (2.13798) | 70,388 | 142 s |
  | T1 dD 0 | 0.8873 (2.55668) | 23,044 | 100 s |
  | T2 dD 0, Q from 0 by degree | 1.0004 (2.88253) | 240,830 | 330 s |
  | T3 dD 1, Q from 0 by degree | 0.7420 (2.13798) | 70,388 | 211 s |

  - **Soundness caveat.** T2 certifies tau = 2.88253, 0.04% above the margin. The margin was checked
    on the first mode, lambda^2 + 2 lambda + 0.1 + e^{-lambda tau} = 0: omega = 0.4947,
    tau* = 2.8814; the second mode crosses at 3.498. Both gates pass the near-margin steps:
    - feasratio 0.867, eta 1.5e-8, relres 8.6e-4;
    - feasratio 0.810, eta 9.4e-9, relres 5.4e-4.

    So a marginally unstable system passes both gates there. Only the relative residual relres
    (about 1e-3, against 1e-6 to 1e-9 on clean certificates) flags it. Reaches within about 0.1% of
    the margin are at the numerical floor of these gates, and reach comparisons are read to 1e-2
    only.
  - On FO2s dD 0 reaches further than the rule's lift. The two families are not nested: a lower lift
    gets a higher weight degree from the reach rule.
- **S23 (07:47-08:43, B1c1 + B1c1g).** B1.1, lpi.ineq.dD = 0 on the 30 two-variable rows, against
  B1c0 (B1_c1_dD0.md). The gain rows 45, 46, 47, 71 were rerun on the committed gain executive
  (B1c1g, D121).
  - Gained by both gates: PDE 61 (gate-refused at light, certified at 197,415 variables).
  - Lost by both gates: PDE 43 (S18, certified -> infeasible).
  - Gates split differently: 38 (now certified by both), 48, 63 (harness unknown, executive optimal).
  - Rows that now fit the gate: 50 (1,585,445 variables; inaccurate, dual unknown); 69 and 70
    infeasible (both documented as not exponentially stable).
  - Unchanged: the gain rows 45-47, 71 stay unknown in 'Q' at 92,587-94,307 variables; 72 goes
    unknown -> infeasible.
  - Totals: variables 7.34M -> 7.83M (1.07x, because rows that now fit are counted); seconds 2,092 ->
    3,344.
  - Verdict: not a tier change. It loses S18; the tightness runs show losses on KISS, S18, S17 and a
    gain on FO2s.
- **09:31.** B1c2 stopped on PDE 65 (the Turing pair). Under dD 0 with the degree raise of Q its
  session held 46.1 GB, the machine had 5.7 GB free, and a further raise would have paged. PDE 65 is
  skipped for the B1 candidates (opts.skip) and recorded as not run for memory. B1c2 resumes from the
  next row.
- **S24 (06:44-09:32, B1 tightness, all seven systems; B1_tightness.md).** Reach as a fraction of
  the known margin, by the harness bisection (tolerance 1e-3 of the margin; read to 1e-2 near a
  margin, D125).

  | system | T0 light | T1 dD 0 | T2 dD 0, Q from 0 by degree | T3 dD 1, Q from 0 by degree |
  |---|---|---|---|---|
  | KISS (34) | 0.9998 | 0.9951 | 0.9998 | 0.9998 |
  | parabolic (37) | 0.9998 | 0.9998 | 0.9998 | 0.9998 |
  | AGE (66) | 0.9881 | 0.9869 | 0.9951 | 0.9951 |
  | S18 (43) | 0.9893 | below 0.5 | 0.9400 | 0.9916 |
  | S17 (72) | 0.678 | below 0.5 | below 0.5 | 0.727 |
  | FO2s (57) | 0.742 | 0.887 | 1.0004 (false, D123) | 0.742 |
  | TL (64) | 0.9869 | 0.9992 | 0.8557 | 0.9869 |

  - **T3** never reaches less than light: higher on AGE, S18, S17; equal on the other four.
    Variables at the reach are 0.40x (KISS) to 2.72x (S17) of light; bisection time 0.56x (para) to
    2.4x (TL, where the degree raises at the failing steps certify nothing more).
  - **T1 (dD 0)** fails S18 and S17 at half the margin, loses 0.5% on KISS, and gains on FO2s (+0.145)
    and TL (+0.012). Not a tier change; the maintainer's tightness concern is confirmed.
  - **T2** loses on S18 and TL and fails S17.
  - **Candidate for the maintainer:** the T3 setting as the light 2-D entry of the direct form:
    sop.deg_Q{2} = 0, sop.raise_Q = 'deg', sop.max_raise_Q(2) = 2. The library run B1c3 gives its
    verdicts at the default parameters.
- **S25 (08:44-09:34, B1c2).** B1.2 on the 29 two-variable rows (PDE 65 skipped): dD 0, deg_Q{2}
  0, raise_Q 'deg', max_raise_Q(2) 2, against B1c0 (B1_c2.md).
  - Gained by both gates: PDE 58 (FO2r, light infeasible -> certified at 225,705 variables) and
    PDE 61 (gate-refused -> certified at 197,415).
  - Lost: PDE 64 (TL). Its passes at Q degree 0 and 1 fail, and the raise to degree 2 builds
    2,089,764 variables, above the size gate.
  - Gate splits: 38 (now certified by both gates, at 1,776,564 variables), 48, 63, 72 (executive
    optimal).
  - Size gate after the raises: 69 (4,456,629), 70 (2,015,785).
  - S18 (43) stays certified at the library parameter (196,513).
  - Totals: variables 1.19x (larger raised programs), seconds 0.85x.
- **S26 (09:32-09:47, study C).** The extreme tier on the 124 one-variable rows (C_extreme_*.md).

  | candidate | certified or converged | infeasible | unknown | lo_fail / hi_pass | variables | seconds |
  |---|---|---|---|---|---|---|
  | C1 extreme as is | 28 | 39 | 34 | 23 | 470,397 | 142 |
  | C2 extreme with the stripped 1-D entries (storage row, Q and free operator [2 2 2]) | 100 | 4 | 16 | 4 | 882,751 | 385 |
  | C3 C2 with dD 0 | 95 | 7 | 18 | 4 | 684,275 | 351 |
  | light (L1ctl) | 105 | 4 | 11 | 4 | 1,462,393 | 654 |

  - C2 and C3 gain 53 and 48 rows over the present extreme by both gates and lose none.
  - Against light, C2 loses 4 certificates (the dual of DDE01, the dual of PDE 21 H-infinity, PDE 30,
    the dual of PDE 7) and the five delay thresholds stripped loses (DDE03 -2.6% to DDE08 -5.5%), at
    0.60x the variables.
  - Against light, C3 loses 9 (those, plus PDE 24, 29, 31, 51, 52 from dD 0) at 0.47x the variables
    and 0.54x the time.
  - For the maintainer: the present extreme certifies 23% of the one-variable rows. A minimal tier
    built from the stripped one-variable entries with dD 0 (C3) certifies 77% at 0.47x the variables
    of light. The two-variable rows of C3 are not run.
- **S27 (09:33-10:28, B1c3).** B1.3, the T3 setting, on the 29 two-variable rows (PDE 65 skipped),
  against B1c0 (B1_c3.md). The setting: dD 1, deg_Q{2} 0, raise_Q 'deg', max_raise_Q(2) 2.
  - Gained by both gates: PDE 58 (FO2r, infeasible -> certified at 302,957 variables) and PDE 61
    (gate-refused -> certified at 841,042).
  - Lost: none.
  - PDE 70 fits the gate and is infeasible (857,016), as documented. PDE 72 (S17) stays a gate split
    (executive optimal at 267,968).
  - Totals: variables 7.34M -> 6.91M (0.94x), seconds 2,075 -> 3,264 (1.57x, mostly the rows that
    now fit: 61 at 365 s, 70 at 485 s).
- **S28 (10:35-10:45, B5, no runs).** The triangular cap per example (B5_capinc_per_example.md): the
  capinc rows of the Fable session's ab_slack runs (D102) against L1ctl, light on the same tree.
  - 64 of 124 rows keep the verdict and the value (1e-4 relative) at fewer variables: candidates
    for st lines in the examples. At ratio <= 0.55: DDE05, DDE06, DDE08, DDE17, NDS07, NDS09, DDF01,
    PDE 1, 9, 13, 17, 68.
  - Gained 5, lost 8 (the rows of D102), value lost beyond 1e-4 on 4 (PDE 20, DDE02, NDS08, DDE11
    dual). Variables 1,462,393 -> 1,112,965 (0.76x).
  - Correction to D102: its base for the 1-D PDE rows predated deg_Q{1} [4 2 2]; its -17% and its
    rows with more variables (PDE 2, 3, 4, 10, 14, 19, 67) mixed the two changes. D102 stands.
- **S29 (10:28-10:45, E3vh).** Veryheavy on PDE 56-72 (results_E3vh.csv, seeded with the lighter
  tiers for the 2-D time projection).
  - Certified: 57, 58, 59, 60, 62, 66, 67, 68, each certified by heavy too; no row heavy misses.
  - Skipped by the projection: 61, 63, 64, 71 (1,268-2,990 s per solve projected). Size gate: 56
    (32.9M variables), 65 (8.08M), 69 (m 72,333), 70 (5.93M). PDE 72 (S17): gate split (fr 1.78,
    executive optimal), as at light.
  - The pass filled (study E): results_all_1118b.csv in the study harness (fill_all.pl), 786 rows;
    complete except heavy on 53-56. Certified or converged: extreme 42, stripped 114, light 118 (of
    158), heavy 128 (of 154), veryheavy 122 (of 158). Every error is a size-gate refusal except
    veryheavy PDE 55 (memory).
- **S30 (09:47-10:58, B4ctl).** Veryheavy on the ten 1-D bisection rows, control of B4. Every value
  and every variable count equals the 11/17 pass (results_all_1118.csv); seconds 1.0-1.6x with the
  second slot busy.
  - Against heavy [6 2 2] (B3q622): equal thresholds on 8 of 10 rows (DDE02-07, NDS05, NDS08) at
    0.15-0.17x the variables of veryheavy (NDS05 145,044 vs 996,660); short on DDE01 (2.0320 vs
    2.0380) and DDE08 (3.179895 vs 3.8956787).
  - Against heavy as is (H1ctl): higher on DDE01, DDE03, DDE08, NDS05, NDS08 at 10.7-42x the variables.
- **S31 (10:45-11:17, C2ext2).** Extreme with the stripped 1-D entries (C2, the best 1-D candidate)
  on the 30 two-variable rows; plan C, step 2. Against the 11/17 extreme rows (control C1ext2 on the
  same tree queued): 14 certified in both, no verdict changed by the harness gate. Variables change on
  PDE 50 only (149,735 -> 305,711, primal and dual H2 rows, both unknown); every other row builds the same
  program; total 6.09M -> 6.40M (1.05x), seconds 1,504 -> 1,883.
- **S32 (11:17-12:20, H2ctl).** Heavy on the 29 two-variable rows (PDE 65 skipped as in B1), the
  file seeded with the light rows of B1c0 for the 2-D time projection: control of H2T3 (the B1
  follow-up, T3 at heavy). All 29 rows equal the 11/17 heavy rows in verdict, value and variables;
  18 certified, PDE 36 and 39 skipped by the projection (5,224 and 3,086 s per solve projected);
  12.72M variables, 3,785 s.
- **S33 (10:58-12:29, B4q866).** Veryheavy with deg_Q{1} [8 6 6] on the ten 1-D bisection rows,
  against B4ctl. Every threshold equal or higher.
  - Higher: DDE01 2.0380078 -> 2.04 (3.01x the variables; the analytic margin, from the factor
    lambda^2 - 0.2 lambda + 1 - e^{-lambda tau}, is tau* = (pi - arcsin 0.28)/1.4 = 2.0412847, and
    the rightmost root at tau = 2.04 has real part -5.1e-4: sound, 0.99937 of the margin); DDE08
    3.8956787 -> 3.9144165 (0.46x the variables; documented 3.9-3.95).
  - Equal on the other eight at 1.13-1.20x the variables. Totals 2,469,828 -> 2,990,286 variables
    (1.21x), 4,246 -> 5,418 s.
  - Not run: the other 1-D rows at veryheavy [8 6 6] (veryheavy 1-D costs about 3.2 h per pass).
- **S34 (12:29-12:36, D1Ph).** Heavy in form 'P' ('Pd' duals) on the 66 DDE and NDSDDF rows, against
  H1ctl (heavy 'Q').
  - No verdict changed by both gates. NDS09 dual, a gate split under 'Q' (harness unknown, executive
    optimal), is certified by both under 'P'.
  - Thresholds: DDE01 2.0180859 -> 2.04, DDE03 1.7161724 -> 1.71785, DDE08 1.42229 -> 3.5471558,
    NDS05 1.5784509 -> 1.5797836, NDS08 0.027388574 -> 0.027523096; DDE02 2.0839844 -> 2.0820312
    (one bisection step, -0.094%).
  - Variables 1,267,915 -> 1,025,620 (0.81x), seconds 554 -> 383.
  - Against light 'P' (D1P): heavy 'P' adds DDE08 (1.42229 -> 3.5471558) and is one step lower on
    DDE02. Light 'P' reaches DDE01 2.04 at 4,428 variables, where veryheavy [8 6 6] uses 283,844.
- **S35 (12:36-12:38, D1Ppde).** Light in form 'P' ('Pd' duals) on the 39 1-D PDE stability rows,
  against L1ctl (light 'Q').
  - Lost by both gates: PDE 10 (unknown), PDE 13.2 and its dual (infeasible), the dual of PDE 7
    (infeasible). Gained: none.
  - Thresholds: PDE 12 from hi_pass 37.5 to 32.582886 (-13.1%); PDE 5 and 6 one bisection step lower
    (-0.098%).
  - Same verdict on the other 32 rows (25 certified or converged), mostly at 0.12-0.21x the variables;
    the duals of PDE 6 and 11 rise 1.46x. Variables 517,760 -> 122,247 (0.24x), seconds 131 -> 113.
  - Unlike the delay rows (D115, D130), 'P' does not keep every certificate of 'Q' on the PDE rows:
    no default for PIETOOLS_PDE; candidates per example only.
- **S36 (12:38-13:13, GallQh, GallPh).** The 21 H-infinity gain rows at heavy, forms 'Q' and 'P'
  ('Pd' duals), against light (GallQ, GallP, D118).
  - Certified by the harness gate: light 'Q' 14, light 'P' 21, heavy 'Q' 16, heavy 'P' 21. Variables
    1.45M, 0.86M, 5.68M, 3.29M; seconds 200, 79, 1,593, 462.
  - Heavy 'Q': PDE 45, 46, 47 unknown (bounds 50.4, 8.51, 306), PDE 71 inaccurate (3.18); gate splits
    on the DDE11 dual (harness refuses) and DDE12, PDE 20 dual, PDE 22 (harness certifies, executive
    optimal at both tiers).
  - Heavy 'P': the 2-D gains 0.11102274, 0.17110568, 0.20669904 against the exact 0.11101, 0.17110,
    0.20668 (within 1.2e-4 relative; light 'P' 0.88%, 0.37%, 0.90% above), PDE 71 0.045943467, PDE 20
    0.5000001 (exact 0.5); elsewhere within 1e-4 of light 'P'.
  - Light 'P' bounds lie within 1.5e-4 relative of heavy 'Q' or below it on every row heavy 'Q'
    certifies, at 0.15x the variables.
- **S37 (12:20-13:21, H2T3).** The B1 follow-up: T3 (dD 1, deg_Q{2} 0, raise_Q 'deg',
  max_raise_Q(2) 2) at heavy on the 29 two-variable rows, against H2ctl (heavy, deg_Q{2} 2).
  - No verdict changed by both gates. PDE 70 fits the gate (857,016 variables) and is infeasible,
    as documented (heavy: size-gate refusal).
  - Gate split: PDE 38, certified by both gates at heavy (deg_Q{2} 2, 744,164 variables, fr 0.793),
    is refused by the harness under T3 (fr -0.121) and accepted by the executive. The raise stops at
    the first pass the executive accepts (passes 2, Q at degree 1, 245,200 variables): the
    executive reads the solver flags only, so a degree ladder ends on an inaccurate pass that the
    tier degree would have skipped.
  - Variables 12.72M -> 10.29M (0.81x), seconds 3,785 -> 3,600.
- **S38 (13:14-13:51, C1ext2).** Extreme as is on the 30 two-variable rows, the same-tree control of
  C2ext2: all 30 rows equal the 11/17 extreme rows. C2ext2 against it: no verdict changed and no gate
  split; only PDE 50 and its dual change (149,735 -> 305,711 variables, unknown in both). S31 holds.
- **S39 (13:51-13:53, D1Ppdeh).** Heavy in form 'P' on the 39 1-D PDE stability rows, against H1ctl.
  - Lost by both gates: the same four as at light (PDE 10, PDE 13.2 and its dual, the dual of PDE 7).
  - No threshold differs (at light 'P' lowered PDE 12 by 13.1% and PDE 5, 6 by one step).
  - Gate split: PDE 41, unknown by the harness in both forms, executive 'optimal' under 'P' with eta
    2.7e-3 (flags-only acceptance; the eta gate of 1e-7 would refuse it).
  - Variables 574,447 -> 323,716 (0.56x), seconds 150 -> 128.
- **S40 (13:22-, D2Pl).** Light in form 'P' on the two-variable stability rows, against B1c0 (light
  'Q'); 19 rows by 14:05, PDE 69 solving after a raise (1,696,905 variables, m 26,907), then paused
  at the end of the window (PAUSE_sop): 70 and 72 not run.
  - Gained by both gates: PDE 61 (gate-refused under 'Q' -> certified at 323,647).
  - Lost by both gates: PDE 43 (S18), 57 (FO2s), 66 (AGE): the complement of the 11/16-17 gap study
    (the direct form certifies them, the storage form does not).
  - Gate splits: PDE 40 (both gates under 'Q'; harness unknown, executive optimal under 'P'); PDE 36
    and 38 (executive optimal under 'Q', infeasible by both under 'P').
  - Variables 5.06M -> 5.81M (1.15x) over the 19 rows.
  - PDE 69 (14:13): infeasible by both gates under 'P' after one raise (1,696,905 variables, m 26,907,
    1,529 s), as documented (not exponentially stable); under 'Q' the size gate refuses it. D2Pl ended
    at the pause (20 rows); PAUSE_sop removed. Every runner stopped at 14:14.
- **S41 (16:45-16:57, merge check).** ndopvar fast-forwarded to origin 7764389f (96c9961b ->
  7764389f, no local commit ahead): DJ 7d0acf2c (size of copvar/cdopvar returns the component
  totals, the block grid by size(P,'block')), AT 12fd1c3e (examples_FI_library_PIETOOLS), DB
  7764389f (DEMO11_poincare_inequality_2D; DEMO11-14 renumbered 12-15; print_convert_summary). The
  uncommitted edits are unchanged (git diff and status identical before and after).
  - M1L (light, the 124 one-variable rows) against L1ctl: 124 of 124 identical on every column
    except the timings (verdict, result, solves, variables, fr, m, iterations, eta, bisection
    steps); 1,462,393 variables in both; 665 s against 654 s.
  - M2L (light, the 14 two-variable rows of B1c0 under 40 s, PDE 46 gain included) against B1c0:
    14 of 14 identical.
  - st_merge_tests (logs M2T_1, M2U_3): Test_copvar 31/31, Test_cdopvar 23/23 and Test_int2sopvar
    pass; Test_qp_copvar_nd runs (case 3; its two warnings are the expected ones). Two stale tests,
    not merge effects: test_dpvar_ops_sop 572 passed, 3 failed ("numeric summand refused"), since AT
    921aed59 (11/17) made c + X legal (X + cI), whose scalar branch returns before the lines of
    7d0acf2c (code path; not rerun on 96c9961b); test_getsol_sop part (i) 9 of 9 at relative error
    below 6e-16, part (ii) stops in its helper cx_PIE2PDEstability, which passes the tier settings
    to poslpivar_settings_sop (reads lpi.pos.int, retired by the storage row of 11/17); the
    executives build the storage through pos_row and POSCOPVAR.
  - examples_FI_library_PIETOOLS 1-17 (plot=0, 13 s in total): every bound at or above its sharp
    constant and within 7e-5 relative of it (15 examples; QP examples 15-17 included); examples 5
    and 6 (no sharp constant) 1.35103 and 2.16395 against 2 and 4 in the cited paper.
  - DEMO11_poincare_inequality_2D: 0.6366 = (b-a)/pi, 42 s.
