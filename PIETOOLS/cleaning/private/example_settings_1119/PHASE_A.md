# Per-example settings, Phase A (MMP, 11/19/2026)

Maintainer, 11/19/2026: "all examples should have optimized settings files. There is no default."
(DECISIONS D140.) Phase A chooses, for each of the 158 library rows, a form and a setting from the
measurements already taken; it runs nothing. Phase B searches below these choices.

## Method

`phaseA.pl` (this folder) reads the results files of the 11/18 studies
(`testing_codes/private/baselines/2026-11-18_studies/private`) and the example-file settings of
EXAMPLE_SETTINGS.md (Sections 3.3, 6, 7, 11.3). The 11/17 pass (`results_all_1118b.csv`) enters
only for a (row, tier) that no study run measured, since the 1-D light runs of 11/17 predate the
`deg_Q` separation of 11/18 01:24.

For each row:
1. Best result over every certified (form, setting): the largest certified parameter (bisection
   rows), the smallest certified value (gains, norms, controllers, estimators), or the verdict
   (stability at the library parameter). Single solves count only with both gates (harness status
   'certified' and executive status optimal or inaccurate); a bisection is decided step by step
   by the harness gate.
2. Form and setting together: the cheapest candidate that reaches the best, in either form,
   ordered by decision variables, then m, then executive passes, then solves, then seconds.
   'P' wins a tie and any case within 10% of the cheapest decision-variable count (D141,
   maintainer 11/19: "P form should be used when it is cheaper and just as good for accuracy.
   Not blindly!"; the first version of this file preferred 'P' wherever it reached the best).
   "Reaches" means within 1e-3 relative, the bisection resolution. On an exact tie otherwise the
   example-file st is kept.

Output: `phaseA_choices.csv` (one line per row) and `phaseA_choices.md` (the table).

## Result

| | rows |
|---|---|
| form 'P' / 'Pd' | 67 (18 of them within 10% of a cheaper 'Q') |
| form 'Q' / 'Qd' | 44 |
| executives without a P/Q form (controllers, estimators, H2) | 26 |
| nothing certified under both gates | 21 |

The 'Q' choices on one-variable rows are often extreme (C1) or extreme with the stripped 1-D entries
(C2, C3), whose 'Q' programs are the smallest measured there (PDE 1: 206 variables against 280 in
light 'P'). No setting below light was measured in 'P', and stripped 'Q' only on the 11/17 tree;
Phase B measures both.

## Consequences that need a decision or a run

1. **The tuned 2-D examples keep 'Q' (D141).** 'P' certifies PDE 33, 35, 37 and 61 at light, but at
   26,545, 26,545, 188,545 and 323,647 variables against 409, 8,009, 16,656 and 197,415 for their
   'Q' settings, so 'Q' stays. The first version of this file moved them to 'P' under the reading
   "'P' wherever it certifies", which the maintainer rejected on 11/19.
2. **DDE08: 'P' certifies less than 'Q'.** 'P' reaches 3.5471558 (heavy); 'Q' with deg_Q{1} [8 6 6]
   reaches 3.9144165, against a documented tauNmax of 3.9-3.95. 'P' above heavy has not been measured.
3. **Seven choices rest on the 11/17 tree** (DDE13, DDE18, PDE 13.2/dual, 31, 32, 49, 49/dual); they
   need a rerun on the current tree.
4. **Fourteen rows take the composite C3 setting** (DDE14, 15, 16, 17/Hinf_control, 20-24,
   25/H2 dual, 26, 26/H2 dual; PDE 10, 13.2): extreme with the stripped 1-D entries and dD 0. In one
   variable only the 1-D entries and the slack knobs act, so a simpler equivalent (stripped with the
   extreme slack knobs and dD 0) is the Phase B check.
5. **21 rows certify under no measured setting:**
   - documented not stable: PDE 69 (not exponentially stable), 70 (not uniformly stable);
   - defects of the library block: DDE19 (legacy output fields, 10/28), DDF03 (no output declared);
   - documented stable or with a documented value, not certified: PDE 4.5, 8 and 8/dual
     (lam < pi^2), 15, 44 (mu < 0.5), 50 and 50/dual (H2, custom settings 0.8502 / 0.8398),
     56 (tau < 1.0347), 65 (0.9 d*), 72 (0.9 r2*; executive gate only on 11/18);
   - no documented value: PDE 16, 36, 39, 41, 53, 54, 55.

## Phase B (started 11/19 01:37)

Cheaper candidates in both forms; the selection of Phase A is rerun over everything measured.
Tags of `st_run` (2026-11-18_studies harness), launcher `st_qrun3.sh`, logs in `logs/` here:
- one-variable rows (124): stripped 'P' (XSP), stripped 'P' with dD 0 (XSPd0), light 'P' with
  dD 0 (XLPd0), stripped as is (XSQ), stripped as is with dD 0 (XSQd0);
- two-variable rows with a Phase A choice: the same five (X2SP, X2SPd0, X2LPd0, X2SQ, X2SQd0);
- PDE 45-47 at heavy 'P' with dD 0 (X2HPd0); DDE08 at veryheavy 'P' (XVP8); the five rows of the
  11/17 tree at veryheavy on this tree (XVH).

Unit costs, measured: light over the 124 one-variable rows 665 s (M1L), stripped 758 s (11/17);
light 'P' over the two-variable 'P' rows 426 s (D2Pl); heavy 'P' on 45-47 90-137 s each (GallPh).
Estimate: about 110 minutes of solver time, about 55 minutes in two MATLAB processes, sessions of
15 minutes (ST_STOPMIN). The three 3-D rows need Sol and keep their measured settings.

## Wiring (not started)

- `examples_PDE_library_PIETOOLS`: take `[PDE_t, st]` from each example file, and pass `st` and the
  form to `PIETOOLS_auto_execute_sop` (as `settings`, `stability_form`), so that no library example
  reaches the tier prompt.
- `PIETOOLS_DDF.m`: read a block `st` as `PIETOOLS_DDE.m` does.
- The tier prompt of `PIETOOLS_auto_execute_sop` stays for PIEs that do not come from a library, and
  `PIETOOLS_PDE.m` keeps 'heavy' for a user system, unless the maintainer extends D140 to them.
