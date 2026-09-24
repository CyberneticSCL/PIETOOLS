# Low-rank 2-D PIE stability certification (proof of concept)

**What this demonstrates.** For 2-D PIE stability, the SDP cone side can be
removed from the scaling path: certifying a Lyapunov operator costs seconds
once a low-rank *face* of the feasible set is known, while assembling the LPI
costs minutes — i.e. **solving time is (much) less than setup time**. On the
demo's heat-equation ladder the same n=1 program that the stock route solves in
1707.3 s certifies here in seconds, and the n=4 member — memory-infeasible for
the stock route (66.8 GiB dense normal matrix on a 63.7 GiB machine) —
certifies in minutes, all of them assembly.

**Read next:** [theory.pdf](theory.pdf) — why certificates are low rank, why
the restricted search is sound, how the face is found (and the alternatives),
and where the method stops working. [GUIDE.md](GUIDE.md) — the full workflow
for certifying *your own* system from the command line.

> **A wider interface and a benchmark suite exist as of 09/23/2026 — see
> [BENCH.md](BENCH.md).** `pielr_solve(PIE,lpi,opts)` mirrors `lpiscript` and
> covers `stability`, `stability-dual`, `l2gain` and `l2gain-dual` across 1-D
> PDEs, 2-D PDEs and DDEs; `pielr_bench` scores each case against an
> interior-point reference solved in the same run and judged by the same gate.
> Two things on this page are affected. First, the acceptance residual there
> is divided by `max|Pop|`, not by `max|Qop|` as item 5 below describes —
> `tests_1d/opcheck.m` records a measured case where the `max|Qop|`
> denominator collapsed and manufactured a threshold that does not exist, so
> the two halves of this package had been accepting on different quantities.
> Second, the ranks quoted in item 2 were measured at `lmit = 400`, and the
> reported rank is **budget-dependent**: on the gain cases it falls by one per
> block between `lmit = 400` and `8000`. `pielr_certify` itself is unchanged
> and remains the baseline the suite measures against.

**The mechanism in 5 lines.**
1. Build the direct-form 2-D stability LPI (the shipped 2-D executive's own
   program, psatz 0): find Gram blocks `X_i >= 0` satisfying equality rows.
2. Measured on this family, certificates exist with per-block rank ~2 (and
   exactly `2n` for `n` replicated states), independent of the operator-basis
   count.
3. Restrict each block to a face `X_i = V_i S_i V_i'` with `V_i` known: the SDP
   collapses to ~20–50 unknowns (seconds), and `S_i >= 0` implies `X_i >= 0`,
   so the restriction can only *lose* feasibility — never fake a certificate.
4. On the replicated state ladder the n=1 face **tensors** (per monomial
   group), so larger `n` needs no discovery and, with the face coefficients
   carried along, no SDP at all.
5. Every certificate is verified at the **operator level**: residual over all
   36 `opvar2d` cells < 1e-6 AND every Gram block PSD. Never an equality
   residual, never a solver flag (both were measured to pass on
   non-certificates).

## The API

```matlab
cert = pielr_certify(PIE)         % PIE from convert(), or a bare struct
cert = pielr_certify(PIE,opts)
```

`PIE` is what `convert(...,'pie')` returns. You can also hand in **raw
operators** without calling convert:

```matlab
PIE = struct();  PIE.T = Top;  PIE.A = Aop;    % your own opvar2d objects
cert = pielr_certify(PIE);                      % vars/dom default from Top
```

Options (all optional): `.settings` (2-D settings struct, default
`set2d_deg(4,[])`), `.rank` (per-block rank to try first), `.face` (skip
discovery: a face from `pielr_tensor`, `face_n1.mat`, or an earlier `cert`),
`.route` (`'bm'` | `'mintrace'` | `'auto'`), `.maxrank` (search ceiling,
default 6), `.gate` (residual gate, default 1e-6), `.refine`, `.seeds`,
`.verbose`.

Returns `cert.ok`, per-block ranks `cert.r`, `cert.op_rel`, `cert.mineig`, the
face `cert.face`/`cert.S`, sizes `cert.Ns`/`cert.unknowns_full`/
`cert.unknowns_face`, timings `cert.t_setup` / `cert.t_discover` /
`cert.t_certify` (kept separate deliberately), `cert.route`, `cert.part`,
`cert.notes` — and prints a one-screen report.

Minimal example (a system that is *not* the demo's):

```matlab
pielr_path();                                   % once per session
pvar s1 s2
x = pde_var('state',1,[s1;s2],[0,1;0,1]);
sys = [diff(x,'t')==diff(x,s1,2)+0.5*diff(x,s2,2)+2*x;
       subs(x,s1,0)==0; subs(x,s1,1)==0; subs(x,s2,0)==0; subs(x,s2,1)==0];
cert = pielr_certify(convert(sys,'pie'));       % cold: discovery + certify
```

A second entry point, `pielr_certify_pos(Ptgt)`, certifies bare positivity of
a GIVEN operator (see GUIDE.md §8 and `demo_poincare_2d.m`).

State ladders: `faceN = pielr_tensor(cert1,n)` replicates a certified n=1 face
to the n-state member of a replicated family; certify with
`pielr_certify(PIEn, faceN)` (no discovery, and with `S` carried no SDP at
all). A previous certificate chains the same way: `pielr_certify(PIE2, cert)`.

## Files

| file | purpose |
|---|---|
| `pielr_path.m` | path setup + shadowing guard; **run first**, edit `PIETOOLS_ROOT` if needed |
| `pielr_certify.m` | the API (above) |
| `pielr_tensor.m` | replicate a certified n=1 face across states, per monomial group |
| `demo_lowrank_2d.m` | the heat-equation showcase (table below) |
| `demo_your_system.m` | commented template: your own PDE, and the raw-T/A entry point |
| `regen_face.m` | regenerates `face_n1.mat` from scratch (~30 min; the .mat is not magic) |
| `face_n1.mat` | shipped certified n=1 heat face (V, S in original units, partition, provenance) |
| [`theory.pdf`](theory.pdf) (`.tex`) | the theory note: why certificates are low rank, why the restricted search is sound, where it stops working |
| [`GUIDE.md`](GUIDE.md) | the usage guide: full command-line workflow for your own system, knobs, failure triage |
| `pielr_certify_pos.m` | bare positivity of a GIVEN opvar2d (one Gram block, every cell matched); same options/report |
| `demo_poincare_2d.m` | worked example for it: the 2-D Poincaré inequality — genuinely degree-hungry |
| `private/` | the measured pipeline: LPI builder, SDP data, face solver, operator-level verifier, BM helpers |

## Running the demo

```matlab
demo_lowrank_2d          % from this folder; ~10 min at the default nlist=[1 2 3]
```

Fresh-run table from this machine (Windows, R2025b, SeDuMi), against the
shipped (regenerated) face:

```
  n | Gram unknowns | setup_s  | solve_s  | solve/setup | op_rel     | psd | route
  1 |        277176 |      5.5 |      1.6 |      0.3021 |  9.291e-07 |   1 | face-replicate
  2 |       1107952 |     24.9 |      6.0 |      0.2395 |  9.291e-07 |   1 | face-replicate
  3 |       2492328 |     65.5 |     14.9 |      0.2273 |  9.291e-07 |   1 | face-replicate
  4 |       4430304 |    278.1 |     47.2 |      0.1697 |  9.291e-07 |   1 | face-replicate
```

(The n=4 row was run through the same API path with `nlist` extended; it is
exactly what `demo_lowrank_2d` produces when you add 4 to `nlist`.)

Reference: the stock full solve of the same n=1 program measured **1707.3 s**
(op rel 1.632e-07); at n=4 the stock route needs a 66.8 GiB dense normal
matrix and is infeasible on a 63.7 GiB machine, while the face route certifies
n=4–6 with no SDP at all (measured op rel 9.291e-07 at every n, the n=1
residual).

## Requirements

- PIETOOLS on the path (`pielr_path` finds it via `pietools_path_update`, or
  set `PIETOOLS_ROOT` at the top of `pielr_path.m`),
- SeDuMi on the path,
- RAM: n=1–3 fits in a few GiB; LPI **assembly** needs roughly 8 / 10 / 13 GiB
  at n = 4 / 5 / 6, where n is the number of replicated PDE states (the demo's
  `nlist`; the certification side stays tiny at every size).

## Regeneration

`regen_face` rebuilds `face_n1.mat` with the package's own discovery
(Burer–Monteiro route) and overwrites the .mat only if the result certifies —
the shipped face is exactly such an output. The deterministic alternative
(`route='mintrace'`: one full SeDuMi solve + eigenvector truncation) was
measured **not** to regenerate this face: SeDuMi exits with numerr=1 and the
trace objective crushes the Lyapunov block to noise, so its truncated faces
fail at every rank <= 6. That is a solution-quality statement about that
route, not a rank floor.

## Caveats (read before quoting numbers)

- **Where "solve < setup" strictly holds:** certification *given a face* —
  always, measured 7–19 s against 10–560 s of setup — and across the state
  ladder by tensoring. For **first contact with a new operator**, discovery is
  a one-time cost (minutes; still several times cheaper than the full solve,
  which took 28 min at n=1 and is memory-infeasible at n=4). The demo table
  and the API report keep discovery and certification as separate columns so
  nobody is misled.
- **Ranks are upper bounds** exhibited by verified certificates at a 1e-6
  operator-residual gate; certificates are approximate in that sense,
  consistent with solver practice. A failed rank is "not reached", never a
  proven floor.
- **Scope of the measurements:** the direct-form 2-D stability LPI at psatz 0,
  on the reaction–diffusion class plus two PIETOOLS library examples
  (Heat_Eq_with_ODE, Reaction_Diffusion, both r* = [2 2]). The r* = 2n law and
  the tensor property are measured regularities of the *replicated* family,
  not theorems for all 2-D PIEs — which is exactly why the API exists: so you
  can probe other cases. First datapoint from the template's anisotropic
  variant: r = [3 3] certified, rank 2 not reached (three BM seeds at
  1.4e-06–4.9e-06, just above the gate). **Retested (CC, 09/22/2026): MARGINAL, cite with
  care.** Those seeds sit only 1.4–4.9x above the gate and were run under
  `bm_lm2`'s then-shipped 400-iteration hard cut, so the budget was re-swept
  with the rank pinned at [2 2] (seeds 11/22/33, lmit 400/2000/8000):

  | seed | 400 | 2000 | 8000 | |
  |---|---|---|---|---|
  | 11 | 2.070e-06 | 1.166e-06 | 1.166e-06 | stagnation exit |
  | 22 | 1.355e-06 | 1.307e-06 | 1.307e-06 | stagnation exit |
  | 33 | 4.888e-06 | 1.082e-06 | **1.000e-06** | still descending |

  The lmit-400 rung reproduces the original 1.4e-06–4.9e-06 range, so this is
  the same measurement. Two seeds hit the windowed stagnation test and stop;
  the third keeps descending and lands at 1.000e-06 against a 1e-6 gate.
  "Rank 2 not reached" therefore stands, but with essentially NO margin, so
  the "+1 per block under anisotropy" departure from r* = 2n is a knife-edge
  threshold effect rather than a robust property of the operator. Do not cite
  it as a rank law.
- **When low rank fails, measured.** (i) If the operator *forces* a
  non-constant multiplier (the delta cell `R22{1,1}` must equal a given
  `a(s1,s2)`, as in weighted integral inequalities or weighted Lyapunov
  functions for variable-coefficient PDEs), the multiplier block's rank is
  exactly the **SOS length** of `a` — measured 1, 2, 3 on a designed ladder,
  with `a = 1+s1²+s2²` forcing the Gram `I₃` — and the whole-Gram floor is
  SOS-length plus ~2 for the integral part. In one spatial variable SOS length
  is at most 2; in two it is unbounded in degree. (ii) A weight that is
  positive on the box but **not globally SOS** (e.g. Motzkin + 0.01) makes the
  LPI infeasible at psatz 0 at *every* `c`, including `c = 0`. (iii) Pure 2-D
  PDE *stability* problems dodge both effects structurally — the negativity
  operator's multiplier cell is zero whenever `T` and `A` carry no multiplier,
  true of all seven PIETOOLS 2-D library examples. (iv) The rank floor rises
  toward full rank as the operating point approaches the LPI's own certifiable
  boundary (measured on the weighted-inequality family: floor 3 at half the
  boundary, 6 at 0.8, none ≤ 6 at 0.95) — low rank is a property of
  certificates with margin, not of borderline ones. **Caveat (CC, 09/22/2026).**
  Those three points were measured with `bm_lm2` at its then-shipped cap of 400
  iterations, which is a hard cut and not a convergence test; see that file's
  header. Raising the budget moved the 1-D reach from 0.41 to 0.999 of `lam*`,
  so the 0.95 entry establishes only that *this search* did not attain rank
  ≤ 6 within that budget — not that no such certificate exists. The
  direction of (iv) is independently consistent with the interior-point rank
  at the boundary (737/744 under a min-trace objective), but the numbers
  3 / 6 / none are not a rank floor and should not be cited as one.
- **Assembly is now the dominant cost.** The point of this package is that the
  cone side has been removed from the scaling path; making `poslpivar_2d` +
  `lpi_eq_2d` assembly cheaper is a separate (open) problem.
