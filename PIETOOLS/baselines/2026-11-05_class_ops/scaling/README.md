# Scaling plots: opvar/dopvar vs sopvar/sdopvar

Data and figures in the style of Figs. 2–4 of D. Jagt, S. Shivakumar, P. Seiler and M. Peet,
"Efficient data structures for representation of polynomial optimization problems: implementation
in SOSTOOLS", IEEE Control Systems Letters 6 (2022) 3493–3498. Measured 11/06/2026 on `ndopvar` at
`8343fd1b` (the measured classes unchanged since `edbd7341`), MATLAB R2025b, MOSEK 11.0.30,
Windows 11, Intel i9-14900KF, 63.7 GB.

| paper | here | file |
|---|---|---|
| Fig. 2: time of s(x;ξ)·p(y) against q | time of the composition D·A of a decision operator D with a constant operator A against q | `fig2_composition_vs_q` |
| Fig. 3: differentiation and substitution against q | addition D + E (disjoint and shared decision-variable lists) and adjoint D' against q | `fig3_ops_vs_q` |
| Fig. 4: parsing time of SOS programs and % of time parsing against problem size | parsing time of the stability LPI and % of time parsing against the degree d of the Lyapunov variable, 1-D and 2-D | `fig4_stability_parse` |

The decision operator classes have no `diff` or `subs`, so Fig. 3 shows addition and the adjoint;
`ops_vs_q.csv` holds every operation of `../README.md` for other panels. Each figure is saved as
`.png` (preview) and `.pdf` (vector).

## Data files

**`ops_vs_q.csv`** – one row per (sweep, q, operation, class): `sweep`, `n`, `d`, `q`, `nnz`,
`op`, `class` (`dopvar` or `sdopvar`), `t` (median seconds), `nrep`, `bytes` (storage of the result).
The sweeps:

| sweep | q | decision structure |
|---|---|---|
| `nnz_fixed` | 13 values, 10³–10⁶ (quarter decades) | each coefficient one decision variable drawn from q; nnz = 110 for every q; the two operands of a binary operation on disjoint lists |
| `nnz_fixed_shared` | as `nnz_fixed` | as `nnz_fixed`, both operands on one list |
| `nnz_4q` | 11 values, 10³–3.2·10⁵ | each coefficient ⌈4q/55⌉ decision variables drawn from q; nnz ≈ 3.9q; disjoint lists |

Operators: n = 1 (scalar operators on L2[0,1]), kernels of degree 4 in s and in t – R0 with 5
monomials, R1 and R2 with 25 each, 55 coefficients – as the paper uses degree 4. Operations:
`plus`, `plus_const` (D + A), `minus`, `uminus`, `scalar` (2.5·D), `mtimes_DA` (D·A), `mtimes_AD`
(A·D), `ctranspose`, `horzcat`, `vertcat`, `blkdiag`, `subsref`.

**`stability_parse_vs_degree.csv`** – one row per (nv, d, series): `nv` (1 or 2 spatial
variables), `d` (Lyapunov degree), `series`, `t_total`, `t_solve` (MOSEK optimizer time), `t_parse`
(= t_total − t_solve; medians over `nrep` calls), `pct_parse` (100·t_parse/t_total), `parse_reps`
(parsing time of every call, first call first), `status`, `ndv` (decision variables of the SDP),
`m` (equality rows), `Ks` (PSD block sizes), `sumKs2` (sum of their squares), `err`.

## Method

**Figs. 2–3.** `../private/bench_class_ops.m`, parts `dvarQ`, `dvarS`, `dvarN` with n = 1, d = 4. Each
operator is built once as a dopvar and converted once (`dopvar2sdopvar`); both classes agree on
every operation to 3.7e-16 (`bench_class_ops('check')`). Time: one warm-up call, then the median of
3–25 calls filling at least 0.3 s.

**Fig. 4.** `../private/bench_lpi_degree.m`. The stability LPI in the direct form, P ≻ 0 and
T'PA + A'PT ≼ 0, for the heat equation x_t = Δx + λx on [0,1]^nv with Dirichlet conditions
(λ = π²/2 in 1-D, π² in 2-D; both stable). Parsing time is the wall time of the executive minus the
SDP solver time reported by MOSEK, as in the paper.

| series | code | classes |
|---|---|---|
| `legacy` | `legacy_PDEstability`, `legacy_PDEstability_2D` (copies of the stock `PIETOOLS_PDEstability(_2D)` at `8b4d1038`, the `pie_struct` check removed) on `pie2opvar(PIE)` | opvar, dopvar (opvar2d, dopvar2d), `lpi_programming` |
| `sop_stock` | `PIETOOLS_stability_sop(PIE,settings,'P')`, `settings.sop.slack = 'stock'`; 1-D only | sopvar, sdopvar, copvar, cdopvar, `lpi_programming_sopvar` |
| `sop` | `PIETOOLS_stability_sop(PIE,settings,'P')`, default slack | as `sop_stock` |

The stock executives no longer run on the current HEAD (the `pie_struct` holds `copvar` operators);
the copies change only the class check, so the LPI and the SDP are those of the stock executives.
Container options: one pass (`max_raise = max_raise_feas = 0`), no witness, no dual read-back.

Degrees, one family per dimension:

- **1-D:** `sosineq_on = 0`, dd1 = dd12 = {d,[d d 2d],[d d 2d]}, dd2 = {d+1,[d+2 d+2 2d+2],…},
  dd3 = {d,[d+1 d+1 2d+1],…}, ddZ = [2d 2d 2d]. At d = 3 this is `lpisettings('veryheavy')`.
  `legacy` and `sop_stock` build the same SDP: equal `ndv` and PSD block sizes at every d, with
  2–13 fewer equality rows in `sop_stock`, since the container layer drops equality rows that are
  multiples of rows already written (`e3a1f667`).
- **2-D:** `lpisettings('light')` with every degree of `settings_2d.LF_deg` multiplied by d (light
  at d = 1); the negativity condition by `lpi_ineq`. The two stacks build different SDPs here (see
  `ndv`, `m`), so the 2-D solve times compare formulations, not data structures.

A call shorter than 10 s is repeated 3 times and the medians are reported. The first points
(q = 10³; d = 1–3 in 1-D) were measured after a discarded warm-up pass: in a fresh session the
first call carries a one-time cost (1-D, d = 1: 0.89 s legacy and 0.23 s container cold, against
0.56 s and 0.05 s warm). At d = 8 the parsing times of three successive calls are 4.04, 4.02 and
4.32 s (legacy) and 0.59, 0.49 and 0.54 s (container): the first call at a new degree is within
20% of the repeats, so the medians also stand for single calls.

## Results

**Fig. 2, composition D·A.**

| sweep | q | dopvar (s) | sdopvar (s) |
|---|---|---|---|
| nnz fixed | 10³ | 0.082 | 3.2e-3 |
| nnz fixed | 10⁶ | 2.95 | 2.9e-3 |
| nnz ≈ 4q | 10³ | 0.113 | 8.5e-3 |
| nnz ≈ 4q | 3.2·10⁵ | 13.98 | 1.90 |

**Fig. 3.**

| operation | sweep | q | dopvar (s) | sdopvar (s) |
|---|---|---|---|---|
| D + E | disjoint lists | 10³ | 3.6e-3 | 9.2e-4 |
| D + E | disjoint lists | 10⁶ | 1.53 | 3.58 |
| D + E | shared list | 10³ | 3.3e-3 | 2.9e-4 |
| D + E | shared list | 10⁶ | 1.75 | 0.051 |
| D' | nnz fixed | 10³–10⁶ | 7.7e-4–8.4e-4 | 3.8e-4–4.5e-4 |

**Fig. 4, stability LPI.**

| nv | d | legacy parse (s) | sop_stock parse (s) | sop parse (s) | legacy % parsing | sop_stock % | sop % | legacy solve (s) | sop solve (s) |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 0.56 | 0.054 | 0.050 | 98 | 82 | 84 | 0.012 | 0.010 |
| 1 | 4 | 1.21 | 0.112 | 0.142 | 90 | 47 | 54 | 0.13 | 0.12 |
| 1 | 8 | 4.16 | 0.52 | 0.52 | 66 | 16 | 28 | 2.10 | 1.34 |
| 1 | 12 | 16.41 | 2.49 | 2.23 | 64 | 19 | 18 | 9.07 | 10.09 |
| 2 | 1 | 3.58 | – | 0.53 | 29 | – | 9 | 8.78 | 5.38 |
| 2 | 2 | 5.41 | – | 2.46 | 12 | – | 6 | 39.4 | 38.5 |
| 2 | 3 | 26.04 | – | 7.70 | 4 | – | 3 | 604 | 222 |

In 1-D the container parses the same SDP in 0.08–0.15 of the legacy time (d = 1–12). In 2-D the
container parses its own SDP in 0.15–0.46 of the legacy time; at d = 3 the solver takes 96–97% of
the total in both stacks.

## Limits

- **Statuses.** Some solves end `unknown` (1-D: d = 9 in `legacy` and `sop_stock`, d = 10 in
  `legacy`, d = 11 in `sop`). Parsing time does not depend on the status; the solve times of those
  rows are times to a stall.
- **2-D.** The 2-D series stops at d = 3, where one legacy solve takes 604 s.
- **Copies.** The legacy series runs copies of the stock executives, changed only in the class check.
- **What the timings cover.** One machine; medians of warm calls. Storage is the `whos` size of
  results, not peak memory.

## Regenerate

From the PIETOOLS root run `pietools_path_update`, then make `../private` the current folder:

    qg = unique(round(10.^(3:0.25:6)));
    Rq2 = bench_class_ops('dvarQ',qg,1,4);  Rs2 = bench_class_ops('dvarS',qg,1,4);
    Rn2 = bench_class_ops('dvarN',qg(qg<=3.2e5),1,4);   save('fig23.mat','Rq2','Rs2','Rn2');
    R1 = bench_lpi_degree(1,1:12,{'legacy','sop_stock','sop'});
    R2 = bench_lpi_degree(2,1:3,{'legacy','sop'});    save('fig4.mat','R1','R2');
    export_scaling_csv('fig23.mat','fig4.mat','../scaling');
    plot_scaling('../scaling','../scaling');

MOSEK must be on the path. The 2-D d = 3 cases take about 14 min; the rest about 25 min.
