# B5: the triangular cap per example (no runs; MMP, 11/18/2026)

Plan B5 (PLAN_1118_studies.md): list the 1-D rows where lpi.ineq.cap_inc 'auto' (the triangular
cap of the Newton polytope note) lowers the variables at the same verdict and value. These are
candidates for the st blocks of examples (D106), not for the tiers.

**Data.** The variant: the capinc rows of the Fable session's ab_slack runs (D102; logs
ab_slack_capinc_{DDE,NDS,PDE1D}.log in that session's scratchpad libpass/logs), light with cap_inc
'auto'. The base: results_L1ctl.csv of the study harness, light on the 124 1-D rows on the same tree
(deg_Q{1} [4 2 2]). Script: testing_codes/private/baselines/2026-11-18_studies/private/b5_list.pl.
A row is a candidate when both verdicts are certified (or both converged), the value is equal
within 1e-4 relative or better in the direction of the row (reach up, bound down), and the variables
drop. Stability rows record the value 1 in the base and NaN in the variant: same verdict.

**Totals (124 rows).** Variables 1,462,393 -> 1,112,965 (0.76x). Candidates 64; gained 5; lost 8;
value lost beyond 1e-4: 4; same verdict and value at equal or more variables: 29; neither
certified: 14.

**Correction to D102.** D102 compared the variant with the 11/17 light rows of results_A.csv for the
1-D PDE rows, built before light moved to deg_Q{1} [4 2 2] (D96). Its variable totals (1,340,899 ->
1,112,965, -17%) and its rows with more variables under the cap (PDE 2, 3, 4, 10, 14, 19, 67) mix the
two changes. Against the same-tree base the cap cuts 24%, and the variables rise only on DDE01, NDS05,
NDS11 and the dual gain of PDE 21 (1.27x, 1.69x, 1.26x, 1.26x). The 8 losses are the same
rows; PDE 68 is certified in both (D102: gained). The decision of D102 (off in every tier) stands.

## Candidates (64), by variable ratio

```
NDSDDF:DDF01             cand       certified   1              -> certified   NaN            ndv    5584 ->    1952 (0.35)
PDE:1                    cand       certified   1              -> certified   NaN            ndv    1396 ->     488 (0.35)
PDE:13                   cand       certified   1              -> certified   NaN            ndv    5584 ->    1952 (0.35)
PDE:17                   cand       certified   1              -> certified   NaN            ndv    5584 ->    1952 (0.35)
PDE:68                   cand       certified   1              -> certified   NaN            ndv   12564 ->    4392 (0.35)
DDE:DDE06                cand       converged   3              -> converged   3              ndv   20349 ->    9949 (0.49)
NDSDDF:NDS07             cand       converged   0.44453613     -> converged   0.44453613     ndv   20349 ->    9949 (0.49)
DDE:DDE05                cand       converged   1.5707         -> converged   1.5707         ndv    5233 ->    2593 (0.50)
DDE:DDE17                cand       converged   0.14833984     -> converged   0.14833984     ndv  132300 ->   65900 (0.50)
NDSDDF:NDS09             cand       converged   0.5994668      -> converged   0.5994668      ndv  130825 ->   64825 (0.50)
DDE:DDE08                cand       converged   1.42229        -> converged   1.42229        ndv    6160 ->    3280 (0.53)
PDE:9                    cand       converged   -0.00024414062 -> converged   -0.00024414062 ndv    2093 ->    1129 (0.54)
PDE:13.2                 cand       certified   1              -> certified   NaN            ndv   10037 ->    6405 (0.64)
DDE:DDE16                cand       certified   0.96293737     -> certified   0.96293702     ndv   86188 ->   59404 (0.69)
DDE:DDE17/Hinf_control   cand       certified   9.3513439      -> certified   9.3513437      ndv   29146 ->   20506 (0.70)
DDE:DDE22                cand       certified   1.3500949      -> certified   1.3500949      ndv   17909 ->   12469 (0.70)
DDE:DDE23                cand       certified   1.8080605      -> certified   1.8080607      ndv   10446 ->    7326 (0.70)
DDE:DDE24                cand       certified   1.0439637      -> certified   1.0440097      ndv   17909 ->   12469 (0.70)
PDE:28                   cand       certified   1.0000533      -> certified   1.0000467      ndv    1180 ->     828 (0.70)
DDE:DDE14                cand       certified   0.39528506     -> certified   0.39528501     ndv    4983 ->    3543 (0.71)
DDE:DDE21                cand       certified   2.3270003      -> certified   2.3270008      ndv    4983 ->    3543 (0.71)
PDE:29                   cand       certified   1.000004       -> certified   1.0000087      ndv    1274 ->     906 (0.71)
DDE:DDE15                cand       certified   0.61043363     -> certified   0.61043365     ndv    5173 ->    3701 (0.72)
PDE:7/stability_dual     cand       certified   1              -> certified   NaN            ndv    1368 ->     988 (0.72)
DDE:DDE20                cand       certified   0.23570228     -> certified   0.23570226     ndv    1626 ->    1210 (0.74)
NDSDDF:NDS03             cand       certified   1              -> certified   NaN            ndv   23540 ->   17716 (0.75)
NDSDDF:NDS12             cand       converged   0.44436562     -> converged   0.44436562     ndv    5885 ->    4429 (0.75)
DDE:DDE04                cand       converged   6.1665521      -> converged   6.1665521      ndv    6196 ->    4708 (0.76)
DDE:DDE07                cand       converged   1.3721821      -> converged   1.3721821      ndv    6196 ->    4708 (0.76)
NDSDDF:NDS01             cand       converged   2.2233267      -> converged   2.2233267      ndv    6196 ->    4708 (0.76)
NDSDDF:NDS02             cand       certified   1              -> certified   NaN            ndv    6196 ->    4708 (0.76)
NDSDDF:NDS10             cand       certified   1              -> certified   NaN            ndv    6196 ->    4708 (0.76)
DDE:DDE03                cand       converged   1.7161724      -> converged   1.7161724      ndv    1712 ->    1324 (0.77)
PDE:10                   cand       certified   1              -> certified   NaN            ndv    2093 ->    1721 (0.82)
PDE:5                    cand       converged   9.8696044      -> converged   9.8696044      ndv    4944 ->    4176 (0.84)
PDE:6                    cand       converged   2.4674011      -> converged   2.4674011      ndv    4944 ->    4176 (0.84)
DDE:DDE26                cand       certified   3.8308551      -> certified   3.8308295      ndv   54996 ->   49044 (0.89)
PDE:49                   cand       certified   0.1016557      -> certified   0.10166147     ndv    1563 ->    1395 (0.89)
DDE:DDE25                cand       certified   3.6746856      -> certified   3.6746586      ndv    1746 ->    1570 (0.90)
DDE:DDE25/H2_norm_dual   cand       certified   3.6723928      -> certified   3.672393       ndv    1746 ->    1570 (0.90)
PDE:21                   cand       certified   0.84187746     -> certified   0.84187751     ndv    3409 ->    3073 (0.90)
PDE:25                   cand       certified   0.81029155     -> certified   0.81029165     ndv    3409 ->    3073 (0.90)
PDE:26                   cand       certified   2.1472936      -> certified   2.1472944      ndv    7893 ->    7125 (0.90)
DDE:DDE09                cand       certified   25.912755      -> certified   25.912633      ndv    1040 ->     948 (0.91)
DDE:DDE09/Hinf_gain_dual cand       certified   25.909445      -> certified   25.909442      ndv    1040 ->     948 (0.91)
DDE:DDE10                cand       certified   0.2364431      -> certified   0.23644312     ndv    3837 ->    3485 (0.91)
DDE:DDE10/Hinf_gain_dual cand       certified   0.23645056     -> certified   0.23645052     ndv    3837 ->    3485 (0.91)
DDE:DDE12/Hinf_gain_dual cand       certified   4.9801891      -> certified   4.9801026      ndv    3837 ->    3485 (0.91)
PDE:23                   cand       certified   23.719575      -> certified   23.717504      ndv     929 ->     841 (0.91)
PDE:24                   cand       certified   8.2137035      -> certified   8.2137124      ndv     929 ->     841 (0.91)
DDE:DDE11                cand       certified   2.8772529      -> certified   2.8774331      ndv    1246 ->    1146 (0.92)
PDE:19                   cand       certified   1              -> certified   NaN            ndv   11875 ->   10975 (0.92)
PDE:14                   cand       certified   1              -> certified   NaN            ndv   22336 ->   20672 (0.93)
PDE:2                    cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:2.2                  cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:2.3                  cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:2.4                  cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:3                    cand       certified   1              -> certified   NaN            ndv   12564 ->   11628 (0.93)
PDE:4                    cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:4.2                  cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:4.3                  cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:4.4                  cand       certified   1              -> certified   NaN            ndv    5584 ->    5168 (0.93)
PDE:67                   cand       certified   1              -> certified   NaN            ndv   34900 ->   32300 (0.93)
PDE:7                    cand       converged   4.6508984      -> converged   4.6508984      ndv    8508 ->    8136 (0.96)
```

## Not candidates

```
PDE:22                   gained     unknown     0.20847412     -> certified   0.20847417     ndv     929 ->     841 (0.91)
PDE:27                   gained     unknown     1.0000223      -> certified   1.0000275      ndv    1180 ->     828 (0.70)
PDE:30                   lost       certified   0.6499404      -> inaccurate  0.70665632     ndv    1394 ->    1094 (0.78)
PDE:31                   lost       certified   0.00025642962  -> inaccurate  0.00020605342  ndv    4387 ->    3331 (0.76)
PDE:51                   lost       certified   1.2699242      -> inaccurate  1.3150868      ndv    2416 ->    1932 (0.80)
PDE:52                   lost       certified   1.0963826      -> inaccurate  1.3164324      ndv    2523 ->    2031 (0.80)
PDE:20                   value-lost certified   0.50081458     -> certified   0.50100408     ndv     929 ->     841 (0.91)
DDE:DDE12                gained     unknown     4.9839982      -> certified   4.9838202      ndv    3837 ->    3485 (0.91)
DDE:DDE13                lost       certified   1.3653282      -> unknown     1.4172884      ndv    4983 ->    3543 (0.71)
DDE:DDE18                lost       certified   0.37791699     -> unknown     0.39343837     ndv   13672 ->   10072 (0.74)
DDE:DDE02                value-lost converged   2.0839844      -> converged   2.0800781      ndv   23540 ->   17716 (0.75)
NDSDDF:NDS08             value-lost converged   0.027388574    -> converged   0.02736167     ndv   24784 ->   18832 (0.76)
PDE:20/Hinf_gain_dual    gained     unknown     0.50000222     -> certified   0.50568005     ndv     929 ->     841 (0.91)
PDE:21/Hinf_gain_dual    lost       certified   0.84188347     -> unknown     0.85393858     ndv    3409 ->    4289 (1.26)
DDE:DDE26/H2_norm_dual   gained     unknown     3.8306734      -> certified   3.8306475      ndv   54996 ->   49044 (0.89)
DDE:DDE01/stability_dual lost       certified   1              -> infeasible  NaN            ndv    5996 ->    4412 (0.74)
DDE:DDE11/Hinf_gain_dual value-lost certified   2.8779832      -> certified   2.8784024      ndv    1246 ->    1146 (0.92)
```

Rows at ratio <= 0.55, where an st line in the example halves the program or better: DDE05, DDE06,
DDE08, DDE17 (bisection and its Hinf control at 0.70), NDS07, NDS09, DDF01, PDE 1, 9, 13, 17, 68.
DDE08 runs 'Q' with deg_Q{1} [8 6 6] in its example block, a configuration not measured with the cap.
