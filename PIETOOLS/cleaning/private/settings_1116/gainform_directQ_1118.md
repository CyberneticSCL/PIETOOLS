# The direct 'Q' form of PIETOOLS_Hinf_gain_sop (MMP, 11/18/2026 07:50)

Request: the Sister session (S16, D118; studies_1118/G_gainform.md) measured the 21 H-infinity gain
rows of the catalog at light in both forms: 'P' certifies 21 of 21 by both gates, 'Q' (the form of
the auto-executor) 14, and the four 2-D rows (PDE 45, 46, 47, 71) are UNKNOWN in 'Q' with numerr 2
in every tier and on HEAD. It asked for the diagnosis of 'Q'. The maintainer is offline (04:40:
"if there is an issue sister agent want you to look at, do it. But I can't make decisions").

## Diagnosis
Read in the executive (not measured): case 'Q' built R from the storage row (tl.storage), Q with
the shape of T at the degrees of R (get_lpivar_degs_sop(R)), and imposed the equality T'Q = R. The
stability executive used the same construction until 11/17 and failed S3 in 2-D with it
(mixed_domain/GAP_RUN_LOG.md); its direct form certifies S3. Inferred cause: in two variables the
storage row is separable and the equality T'Q = R admits no feasible R. Measured support: replacing
the construction certifies the four 2-D rows (below).

## The candidate (OUT OF THE TREE since 07:54, at the request of the Sister session, whose controls run
## on the committed tree: git patch and full copies (.m.txt) in gainform_candidate/ here; shadow copies in
## the Fable scratchpad gainq_direct/ (dwbase 0, the tested one) and gainq_dw/ (dwbase = dw_Q, queued))
- case 'Q': Q = lpivar_cdopvar(prog,T,qdeg), qdeg = sop.deg_Q{nv} ([] = lpi.var.deg); the equality
  T'Q = Q'T; R = (T'Q + Q'T)/2 >= 0 by lpi_ineq_sop with a slack sized from its own support (no
  'like'), dw = dwbase + dw, cap from R; the KYP slack likewise sized from its own support.
- dwbase = 0 in form 'Q' (the choice of the stability executive). Under test (08:00): the variant
  that keeps dwbase = sop.dw_Q (1 at light) as the former 'Q' did; results_gainQdirect_dw1.csv.
- 'Qd', 'P', 'Pd' unchanged. 'Qd' is routed to 'Pd' for nv >= 2 since 10/29/2026, so the dual
  form has no 2-D exposure; in 1-D its rows are identical to the former run.
- Help text, header entry; (was) lines keep the former construction.

## Measured: ab_gainq.m (Fable scratchpad libpass/ab), light, all 21 gain rows, the harness gate
## |OPT_STATUS - 1| <= 0.3; results in testing_codes/private/baselines/2026-11-17_tiers/private/
## results_gainQdirect.csv; log libpass/logs/ab_gainq_new.log; 76 s wall
Against the former 'Q' (results_GallQ.csv) and 'P' (results_GallP.csv) of
testing_codes/private/baselines/2026-11-18_studies/private. Entries: status, gamma, variables.

| row | former 'Q' | direct 'Q' (dwbase 0) | 'P' | exact |
|---|---|---|---|---|
| PDE 20 | certified 0.50081458, 929 | unknown 0.53800611, 1585 (MOSEK optimal, eta 1.1e-10, OPT_STATUS 0.52) | certified 0.50000019, 1177 | 0.5 |
| PDE 20 dual | unknown 0.50000222, 929 (OPT_STATUS 0.657) | identical | certified 0.50000031, 1177 | 0.5 |
| PDE 21 | 0.84187746, 3409 | 0.84200387, 5953 | 0.84187719, 4345 | |
| PDE 21 dual | 0.84188347, 3409 | identical | 0.84187744, 4345 | |
| PDE 22 | unknown 0.20847412, 929 (gate split) | certified 0.20833338, 1585 | 0.20833341, 1177 | |
| PDE 23 | 23.719575, 929 | 23.705799, 3997 | 23.705775, 1177 | |
| PDE 24 | 8.2137035, 929 | 8.2141221, 1585 | 8.2148936, 1177 | |
| PDE 25 | 0.81029155, 3409 | 0.81021217, 5953 | 0.81021921, 4345 | |
| PDE 26 | 2.1472936, 7893 | 2.1456966, 22857 | 2.1457031, 10041 | |
| PDE 45 | unknown, 351307, 52.5 s | certified 0.11320384, 74830, 9.7 s | 0.111982, 200986, 18.4 s | 0.11101 |
| PDE 46 | unknown, 351307, 42.0 s | certified 0.1983716, 74830, 7.3 s | 0.17173673, 200986, 16.7 s | 0.17110 |
| PDE 47 | unknown, 351307, 41.6 s | certified 0.20671187, 74830, 8.6 s | 0.20853507, 200986, 11.9 s | 0.20668 |
| PDE 71 | unknown, 354659, 46.5 s | certified 0.050699177, 75225, 8.2 s | 0.045944081, 203522, 15.1 s | 0.0459 (the LRM record) |
| DDE09 | 25.912755, 1040 | 25.93329, 1746 | 25.907746, 1292 | |
| DDE09 dual | 25.909445, 1040 | identical | 25.907747, 1292 | |
| DDE10 | 0.2364431, 3837 | 0.23644302, 6581 | 0.23644313, 4789 | |
| DDE10 dual | 0.23645056, 3837 | identical | 0.23644327, 4789 | |
| DDE11 | 2.8772529, 1246 | 2.8667247, 2025 | 2.8666078, 1515 | |
| DDE11 dual | 2.8779832, 1246 | identical | 2.8666082, 1515 | |
| DDE12 | unknown 4.9839982, 3837 (gate split) | certified 4.9792396, 6581 | 4.9784482, 4789 | |
| DDE12 dual | 4.9801891, 3837 | identical | 4.9784485, 4789 | |

Reading:
- Certified by the harness gate: direct 'Q' 19 of 21; former 'Q' 14; 'P' 21.
- Gained against the former 'Q': the four 2-D rows at 0.21x the variables and 0.2x the time, and
  the two gate splits DDE12 and PDE 22.
- Lost against the former 'Q': PDE 20 only. Its program solves to optimality (eta 1.1e-10) at a
  bound 7.6% above the exact 0.5; the executive accepts it, the harness gate refuses it at
  OPT_STATUS 0.52. The raise ladder does not act on an optimal solve, so the base raise decides
  the bound: the former 'Q' solved at the base raise dw_Q = 1. The dwbase variant (above) tests
  this.
- Against 'P' in 2-D: the direct 'Q' is tighter on PDE 47 (-0.9%) and looser on 45 (+1.1%), 46
  (+15.5%) and 71 (+10.3%), at 0.37x the variables and 0.5-0.7x the time.
- Against 'P' in 1-D: every certified primal row within 0.1% of 'P'; the direct 'Q' costs 1.3-2.3x
  the variables of 'P' and 1.7-4.3x those of the former 'Q' (own-support slacks and Q at
  deg_Q [4 2 2] against Q at the degrees of R). In 1-D the ordering in variables is
  'P' < former 'Q' < direct 'Q' on every primal row.

## Open (the maintainer's decisions)
- Whether the candidate replaces the former 'Q' (and with which base raise), and whether the
  auto-executor and the entry scripts pass 'P' for gains (D118's question). Nothing committed.

## The base-raise variant (dwbase = sop.dw_Q = 1 kept, as the former 'Q'; shadow gainq_dw/;
## ab_gainq_dw.m; results_gainQdirect_dw1.csv; log libpass/logs/ab_gainq_dw_dw1.log; 207 s wall; 08:04)
Against the dwbase-0 candidate (results_gainQdirect.csv) and 'P'. Entries: status, gamma, variables, seconds.

| row | direct 'Q', dwbase 0 | direct 'Q', dwbase 1 | 'P' | exact |
|---|---|---|---|---|
| PDE 20 | unknown 0.53800611, 1585 (OPT_STATUS 0.52) | certified 0.53800604, 3141 (OPT_STATUS 0.844) | 0.50000019, 1177 | 0.5 |
| PDE 45 | 0.11320384, 74830, 9.7 s | 0.11100942, 538621, 31.9 s | 0.111982, 200986, 18.4 s | 0.11101 |
| PDE 46 | 0.1983716, 74830, 7.3 s | 0.17107346, 538621, 45.7 s | 0.17173673, 200986, 16.7 s | 0.17110 |
| PDE 47 | 0.20671187, 74830, 8.6 s | 0.20668235, 538621, 30.9 s | 0.20853507, 200986, 11.9 s | 0.20668 |
| PDE 71 | 0.050699177, 75225, 8.2 s | 0.045942484, 541973, 58.0 s | 0.045944081, 203522, 15.1 s | 0.0459 (LRM) |
| DDE11 | certified 2.8667247, 2025 | unknown 2.8667167, 3753 (OPT_STATUS 0.553) | 2.8666078, 1515 | |
| other 1-D primal rows | | same bounds to 5e-5 relative at 1.6-2.0x the variables | | |
| dual rows | | identical ('Qd' unchanged) | | |

Reading:
- PDE 20 keeps the bound 0.538 at the base raise 1: the harness gate now accepts it (OPT_STATUS
  0.844) but the 7.6% loss against the former 'Q' (0.50081) stands. The hypothesis "the base raise
  decides the bound" is refuted for this row; the loss lies in the parameterization of Q (the shape
  of T at deg_Q [4 2 2]) against the former Q at the degrees of the storage-row R. Under test:
  deg_Q raised on PDE 20 alone (ab_gainq_p20.m).
- In 2-D the base raise 1 reaches the exact gains to five digits on all four rows, tighter than 'P'
  (45 -0.9%, 46 -0.4%, 47 -0.9%, 71 -0.004%), at 2.7x the variables of 'P' and 7.2x those of the
  dwbase-0 candidate, 2-4x the time of 'P'.
- DDE11 loses the harness gate at the base raise 1 (bound unchanged to 3e-6).
- Three options for 2-D gains at light, from the measurements: the dwbase-0 candidate (cheapest,
  loosest: +1 to +16% against 'P'), 'P' (within 0.4-0.9% of exact), the dwbase-1 candidate (exact
  to five digits, dearest).

## PDE 20 (transport with disturbance, exact gain 0.5): degree against structure (ab_gainq_p20.m,
## 08:44; results_gainQdirect_p20.csv; log libpass/logs/ab_gainq_p20_p20.log; 31 s wall)
Direct form at dwbase 0, light, sop.deg_Q{1} swept, through two shadows: Q with the shape of T (the
candidate) and Q of full shape at the same degrees (the former constructor). The two columns are
identical at every degree (same variables, same bound), so the shape of T does not bind on this PIE.

| deg_Q{1} | status | gamma | variables | m |
|---|---|---|---|---|
| [4 2 2] (light) | unknown (OPT_STATUS 0.52) | 0.53800611 | 1585 | 138 |
| [6 2 2] | unknown (0.36) | 0.53800511 | 2631 | 159 |
| [8 2 2] | unknown (0.48) | 0.53800459 | 4513 | 197 |
| [4 4 4] | certified | 0.50000012 | 6525 | 288 |
| [6 4 4] | certified | 0.50000012 | 6527 | 288 |
| [6 6 6] | certified | 0.5 | 24891 | 530 |
| [8 6 6] | certified | 0.5 | 24893 | 530 |
| former 'Q' (committed) | certified | 0.50081458 | 929 | 124 |

Reading: the 7.6% loss is the degree of the integral kernels of Q (entries 2 and 3 of deg_Q), not
its multiplier degree (entry 1) nor its shape nor the base raise. Kernel degree 4 reaches the exact
gain to 2.4e-7 at 6,525 variables; 'P' at light reaches it to 3.8e-7 at 1,177. In 1-D the direct
'Q' needs kernel degrees above light's [4 2 2] on this row and stays 5.5x dearer than 'P'.
