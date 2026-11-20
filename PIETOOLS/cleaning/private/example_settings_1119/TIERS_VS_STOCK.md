# Present settings tiers against the stock tiers (MMP, 11/19/2026; no new runs)

Script: `tiers_vs_stock.pl`. Stock: stock executives and stock settings, baseline of 10/28/2026
(`testing_codes/private/baselines/2026-10-28_example_library/stock_results.csv`). Present: lpisettings_sop
tiers through the container executives in the library default forms ('Q' for stability and gains), current
tree (11/18-11/19 control passes: light M1L, B1c0; heavy H1ctl, H2ctl; stripped XSQ, X2SQ, E1str; extreme
C1ext, C1ext2; veryheavy B4ctl, E3vh and the 11/17 pass). The 138 rows of the baseline; the harness gate on
both sides (the stock results have no executive gate). Decision variables = cone size nx on both sides.

| tier | certified present / stock | present only | stock only |
|---|---|---|---|
| stripped | 108 / 21 | 87 | 0 |
| light | 112 / 101 | 20 | 9 |
| heavy | 119 / 104 | 16 | 1 |
| veryheavy (mostly 11/17 tree) | 115 / 106 | 11 | 2 |
| extreme | 35 / 14 | 22 | 1 |

Where both certify: light thresholds higher 6, equal 12, lower 1; values lower (better) 10, equal 16, higher
0. Heavy thresholds higher 5, equal 16, lower 1; values better 3, equal 30. The lower threshold is NDS08 in
both tiers: 0.027389 against 0.027523 (-0.49%).

Cost: one variable (120 rows), light median 2.58x the stock variables (total 3.39x), heavy 2.03x (2.72x),
seconds +16%; two variables, light total 0.75x (time 911 -> 231 s), heavy total 0.14x (7 rows; time 2,526 ->
109 s). The one-variable excess is the default form 'Q' with deg_Q{1} [4 2 2]; 'P' at light reaches the same
delay thresholds at 1.4-4.6x fewer variables than any 'Q' setting (EXAMPLE_SETTINGS.md Section 5).

## Complexity at a fixed accuracy (cost_at_fixed_accuracy.pl, .csv)

Tier levels of the two families need not be comparable, so per row a common target: the best accuracy both
the stock tiers and the present tiers reach (stability: a certificate; bisection: the lower of the two best
thresholds; gains and norms: the higher of the two best values; 1e-3 relative). Cost of a family: the fewest
decision variables over its measurements that reach the target (stock: its six tiers; present tiers: the
five tiers in default forms; tuned: every present measurement, including the Phase B variants and the knob
search finals). 117 rows reach a common target; 5 rows are certified by stock only, 5 by the present tiers
only, 11 by neither.

| rows | n | present tiers / stock: median (IQR), geometric mean, rows < 1 | tuned / stock: median (IQR), geometric mean, rows < 1 |
|---|---|---|---|
| 1-D gains, norms | 41 | 1.19 (1.06-1.54), 1.08, 10 | 0.49 (0.25-0.76), 0.40, 38 |
| 1-D stability | 48 | 1.85 (0.97-5.04), 1.90, 20 | 0.26 (0.13-0.42), 0.22, 45 |
| 1-D thresholds | 22 | 2.93 (1.82-4.52), 2.46, 2 | 0.27 (0.11-0.35), 0.24, 19 |
| 2-D stability | 6 | 4.40 (1.54-8.05), 3.32, 1 | 0.51 (0.09-0.71), 0.24, 6 |
| all | 117 | 1.66 (0.97-2.94), 1.68, 33 | 0.33 (0.13-0.61), 0.28, 108 |

Tuned costlier than stock (1.03-1.9x): PDE 5, 6, 11 (thresholds), 24, 25, 29 (gains), 13.2 and its dual,
7 dual. The stock family picks its cheapest of six tiers per row; the tuned family the cheapest of all
measured present settings (up to about 60 per row), so the ratio measures what example-specific settings
buy, not the efficiency of the code at one setting. Time is not compared (one stock run per row on 10/28).

## The controlled comparison of the settings report (11/16/2026, Section "The tiers against the baselines")

Same executive, form 'P', same raise loop and slack terms; only the storage replaced (stock storage
against the tier rows):
- One variable: the light row at the stock lengths is the stock storage: the same variables, rows and
  verdicts on all eight rows (the uniform row without incomplete moments is 23% smaller and loses PDE 3).
- Two variables, light: the uniform and separable rows match the stock sizes within 0.6%, the same
  verdicts, gains within 1e-4 of stock; on the fin2 ladder the stock gain within 2% of its gap at the
  stock size within 0.3%. On the mixed suite the rows certify S3 and S8, which the stock storage leaves
  UNKNOWN after its raise; S18 only under the separable row at d2 = 3 (veryheavy), where stock is infeasible.
- Two variables, heavy: the separable row is stock-heavy-sized (within 0.2%) and certifies S3, but its fin2
  gap is 1.5-1.9e-4 against the stock heavy 3.7e-5 (the stock cap tables admit monomials that lengths and
  per-direction caps do not; a recorded limit of the report).
So the storage rows are stock-equal in one variable and stock-sized in two; the one-variable cost excess
of the table above comes from the default executive form ('Q', deg_Q{1} [4 2 2]) and its slack, not from
the storage.

Stock only at light: PDE 45, 46 (2-D H-infinity gains; the auto-executor runs gains in 'Q', which fails in two
variables, 'P' certifies within 0.4-0.9% of exact), PDE 22, 27, DDE12, the H2 duals of DDE26 and PDE 49,
NDS11 dual, DDE19 (legacy output fields: the stock value is meaningless). At heavy: DDE11 dual.
