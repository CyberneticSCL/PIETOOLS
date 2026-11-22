# Problem cases – settings variations and trades (MMP, 11/20/2026)

Maintainer, 11/20/2026 ~05:00: "create and run a plan for trying different settings variations and trades to
certify these problem cases. Run autonomously for 11 hours." Window: 11/20 05:10 – 16:10.

## Cases

The examples no tested setting certifies (PHASE_C.md §4, §9, §12, §13):

| Case | Signature so far | Not tried before this run |
|---|---|---|
| 41, heat with delayed boundary input (1-D) | unknown on 12 starts, feasibility ratio −0.80 to 0.43; SeDuMi inaccurate (0.32, 0.37) | margins eppos, epneg (1e-4 and 0 in every tier); SDPT3; Q degree above [6 6 6]; storage trades |
| 44, wave with partly delayed boundary damping (1-D) | infeasible with certificate, 8 of 9 starts; shifts to c = 0.05 infeasible | Q degree [8 8 8], [10 6 6], [2 6 6]; margins; 'P' storage above veryheavy; SDPT3 |
| X3, the system of 44 on two variables | infeasible at every tier | Q degree raises in both entries; stock 2-D storage in 'P'; margins |
| 36, 2-D telegraph | MOSEK optimal at ratio 0.375 / 2.39 (53k), unknown at 744k and 1.78M | margins; deg_Q{2} [1 1 1], [2 2 2], [0 1 1]; psatz{2} variants; SDPT3; stock 2-D storage in 'P' |
| 39, 2-D Navier–Stokes | ratio 2.64 at 317k, 7.14 at 649k (MOSEK optimal); 'P' stock light unknown 0.498 | margins; deg_Q{2} [1 1 1], [2 2 1]; psatz{2} 'faces1'; SeDuMi; stock 2-D storage at heavy |
| 69 (Batty), 70 (Beale), 2-D | infeasible with certificate on every start | stock 2-D storage in 'P'; deg_Q{2} raises under the triangular cap; psatz{2} 'faces1'; margins |
| 56, 2-D heat with 3-D delay state | infeasible at every workstation-sized program | 'P' stripped at the default dD; 'Q' stripped with the triangular cap; margins |
| 50 dual, H2 | 'cco' infeasible (96k, 1.66M); 'c' unknown at every tier | 'cco' at the default dD; margins; 'c' heavy with dD 0 |
| 55, 3-D | Sol job 2067511 running since 03:02 | – (collect the result) |

Trades: the grid pairs a raise of one family with a cut of the other – Q degree up with dD 0 or the triangular cap,
storage up with dD 0, slack up (dD 2) with storage or Q degree down.

Margins: every tier sets eppos = 1e-4 and epneg = 0. A larger eppos rescales the storage; epneg > 0 asks for
exponential decay, a stronger claim, which the runs report as such. The inconclusive rows (41, 36, 39) are the
targets; for 69 and 70 (documented not exponentially / not uniformly stable) epneg stays 0.

## Infrastructure

- Clean worktree `C:/wt1120` at aad9606c (the tree every setting was verified on); the main checkout carries the
  edits of the release comment pass. Anything that certifies is rerun on 41c9210c (HEAD of the main checkout since
  11/20 ~04:00: Declan's dpvar subsref and compress) before it is reported.
- Driver `pc_run.m` (harness folder of the worktree): lb_setup first (restoredefaultpath, the worktree alone on
  the path), then diag_amt (one build and solve per job) or ex_search. Job grid `pc_jobs.csv` (generator
  pc_jobs_gen.pl, here); H2 dual starts `pc_starts_50d.csv`. Sanity run `pc_sanity.m`: 13 functions resolve to the
  worktree, the statement split keeps column vectors, one job of 41 ran.
- `ex_search.m` (worktree copy): apply_mods splits a start at a ';' outside brackets only (PHASE_C.md §13);
  TRADE3 reruns the trade pass on the 80 one-variable rows TRADE1 missed (an optimization, not a certification).
- Gate: 2e6 variables, m 45,000; certified = harness gate |fr − 1| <= 0.3, no numerical error, no infeasibility
  flag (lb_exec). Sessions of 60 minutes (ST_STOPMIN), at most two MATLAB processes; the two-dimensional
  waves run exclusive. No Sol job above 10 CHE (standing rule); the uncapped PDE 56 program stays held.

## Queues and prices

| Queue | Tags | Jobs | Price (measured units) |
|---|---|---|---|
| 1 | W41, W44, WX3, EXS_50d, TRADE3 | 42, 17, 9, 6 starts, 80 rows | 41: 1–11 s per job; 44: 18–45 s at 0.09–0.27M, raises unpriced; X3 10–60 s; 50 dual 5–412 s; TRADE1: 10.5 min for 35 rows. About 1.5 h |
| 2 | W36, W39 | 15, 9 | 36: 5 s at 53k, 170 s at 754k; 39: 49 s at 317k, SeDuMi unpriced. About 40 min |
| 3 (exclusive) | W2S, W2Q, W56 | 4, 10, 4 | 69: 49–160 s at 0.34–1.0M; 70: 39–43 s at 0.4–0.5M; 65 stock: 467 s at 1.64M; 56: 178 s at 840k. About 1 h |

Follow-ups, in order: (1) a certified job → ex_search descent from it, the library path with its st in the
worktree copy of the example, the rerun on 41c9210c; (2) directions that move the feasibility ratio toward 1 →
finer grids along them; (3) the 2-D plate (Jagt–Peet; not in the catalog: a direct build); (4) PDE 55 result.
Example files in the main checkout are not edited during the run (the comment pass is editing them): certified
settings are reported with their statements, for the maintainer to apply.

Records: RUN_LOG.md (here), logs/, results copied from the worktree harness at the end.
