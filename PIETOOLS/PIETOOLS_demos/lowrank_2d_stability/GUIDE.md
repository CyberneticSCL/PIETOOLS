# Guide: certifying stability of your own 2-D system

The full workflow, start to finish, defining the PIE from the MATLAB command
line. Why any of this works — and where it stops working — is in
[theory.pdf](theory.pdf); this page is only *how*.

## 0. One-time setup

```matlab
cd <this folder>
pielr_path();        % once per MATLAB session
```

`pielr_path` finds PIETOOLS via `pietools_path_update` (or set
`PIETOOLS_ROOT` at the top of `pielr_path.m`), requires SeDuMi, and errors
loudly if anything on the path shadows the routines this package relies on.
Run the certifying calls **from this folder** — 2-D conversion has been
observed to fail from directories cluttered with loose `.m` files.

## 1. Define the PDE at the command line

Standard PIETOOLS command-line input. Example — anisotropic diffusion with a
reaction term on the unit square, Dirichlet on all four edges:

```matlab
lam = 2;
pvar s1 s2
clear stateNameGenerator            % pde_var keeps counters between definitions
x   = pde_var('state',1,[s1;s2],[0,1;0,1]);
sys = [diff(x,'t')==diff(x,s1,2)+0.5*diff(x,s2,2)+lam*x;
       subs(x,s1,0)==0;  subs(x,s1,1)==0;
       subs(x,s2,0)==0;  subs(x,s2,1)==0];
```

Notes that save time: the second argument of `pde_var` is the number of
states; the domain is a 2×2 matrix, one row per spatial variable; boundary
conditions are `subs(...)==0` rows in the same system list; and the
`clear stateNameGenerator` line matters when you redefine systems in one
session.

## 2. Convert to a PIE

```matlab
PIE = convert(sys,'pie');
```

## 3. Certify (cold call)

```matlab
cert = pielr_certify(PIE);
```

**What is being searched for.** Two positive PI operators, each carried by its
own Gram block, with **separate polynomial degrees**:

- **LF** — the Lyapunov operator `P ⪰ 0` itself; its degree is the *Lyapunov
  degree* of `opts.settings`;
- **eq** — the negativity certificate `D ⪰ 0` in `T'*P*A + A'*P*T = -D`; its
  degree is set separately and higher (the operator products raise degree, so
  the certificate needs the richer basis).

Each block gets its own face and its own rank — that is why the report shows
two of everything: in `r = [3 3] of [10 424]`, the first entry is the LF
block, the second the eq block, and the eq Gram is the large one.

This does four things, timed separately in the report:

| phase | what happens | typical cost |
|---|---|---|
| setup | assemble the stability LPI (same program as `PIETOOLS_stability_2D`, psatz 0) and extract the SDP data | seconds |
| discovery | search for a low-rank face: factored descent from a few seeds, starting at the ESTIMATED rank 2·n₂ (the measured law) and escalating to `opts.maxrank` | **minutes** (measured 11–25 min on two cold problems) |
| certification | solve the ~10–50-unknown restricted SDP on the found face, then shrink it per block | seconds |
| verification | rebuild the PI operators from the candidate and check the operator identity over all 36 parameter cells, plus PSD of every block | seconds |

A one-screen report prints at the end. The line that matters:

```
verdict: CERTIFIES   op rel = 9.897e-07 (gate 1e-06)   PSD = 1
ranks  : r = [3 3] of [10 424]   (Gram unknowns 90155 -> 12, 7513x fewer)
```

**What the input actually is.** There is nothing special about the `convert`
output: `pielr_certify` reads exactly two fields of its first argument —
`PIE.T` and `PIE.A`, the `opvar2d` operators of `d/dt (T x_f) = A x_f` — and
ignores everything else (inputs, outputs, disturbances; the routine certifies
autonomous stability only, the LPI `T'*P*A + A'*P*T ⪯ 0` with `P ⪰ 0`). So
any struct carrying those two operators is a valid input, whether `convert`
produced it or you assembled it yourself:

```matlab
raw = struct();  raw.T = Top;  raw.A = Aop;   % your own opvar2d objects
cert = pielr_certify(raw);
```

The two operators must have matching dimensions on a common domain — exactly
what `convert` returns — and the spatial variables/domain are read off `T`
itself (`T.var1`, `T.var2`, `T.I`; override with `.vars`/`.dom` fields if you
must). The package is 2-D-specific: 1-D `opvar` inputs are refused with a
clear error. `demo_your_system.m` runs both a converted and a hand-assembled
input on the same system.

## 4. Reading the result

- `cert.ok` — a verified certificate exists at `cert.r`. **Ranks are upper
  bounds**: they are exhibited by a certificate, never proven minimal. A rank
  the search failed at is reported as *not reached* — that is a statement
  about the search, not about the problem.
- `cert.op_rel`, `cert.mineig` — the verification numbers; the gate is
  `op_rel < 1e-6` *and* every block PSD. Neither alone suffices (both
  single-condition acceptances were measured to pass non-certificates).
- `cert.t_setup / t_discover / t_certify` — quote these separately.
  "Solving < setup" refers to certification; discovery is a one-time cost.
- `cert.notes` — accumulated fine print for this run; read it.

## 5. Reuse a previous certificate — this is where it pays

A previous certificate **is itself a valid input** to the certify routine:
pass it straight in as the second argument and discovery is skipped entirely —
the new call is certification only, seconds instead of minutes.

```matlab
cert2 = pielr_certify(PIE2, cert);     % seconds, no discovery
```

Two situations where this pays:

- **Re-certifying after a change** — a tweaked parameter, a coefficient, a
  richer degree setting. If the old certificate's face still contains a
  certificate for the new problem, you get it in seconds; if not, you get an
  honest failure (the restriction can only lose feasibility, never fake a
  certificate) and you fall back to a cold call. **Expect the failure to be
  common across parameter changes** — the face encodes the operating point.
  Measured on the demo's heat equation, λ = 2 → 1.9: the shipped face fails
  at op rel 1.7e-2 (an O(1) miss, not a near-thing), while a cold call
  re-certifies λ = 1.9 at the *same* rank [2 2] (op rel 8.2e-07, discovery
  ~3 min). A face failure says nothing about the new system's stability —
  only that this face doesn't contain its certificate.
- **Replicated states** — for a family of `n` copies of a certified 1-state
  system, the face is *constructed*, not searched:

```matlab
faceN = pielr_tensor(cert1, n);        % per-group replication, layout checked
certN = pielr_certify(PIEn, faceN);    % measured: n=1 residual reproduced at n=2..6
```

(`pielr_tensor` also takes a previous cert directly.) Certificates can be
saved and reloaded like any struct — `save mycert.mat cert`, and later
`load mycert.mat; cert2 = pielr_certify(PIE, cert);` — which is exactly how
the demo ships its heat-equation face.

## 6. The knobs (`opts`, all optional)

| field | default | use it when |
|---|---|---|
| `.settings` | `set2d_deg(4,[])` | change either of the two degrees — see below |
| `.rank` | 2·n₂ (estimated) | discovery STARTS here; lower ranks come from the refine shrink |
| `.face` | — | skip discovery (see §5) |
| `.route` | `'auto'` | `'bm'` factored descent / `'mintrace'` deterministic fallback |
| `.maxrank` | 6 | raise if nothing certifies by rank 6 |
| `.gate` | `1e-6` | don't loosen without a reason you can defend |
| `.seeds` | `[11 22 33]` | add seeds when the search near-misses |
| `.refine` | auto | per-block shrink after acceptance |

**About `.settings` — the two degrees are separate, and both can be changed.**
Any 2-D settings struct is accepted, and it carries the two degrees
independently: `LF_deg` (the Lyapunov operator's basis) and `eq_deg` (the
negativity certificate's basis). Two ways to supply one:

- **The convenience helper** `set2d_deg(Dup,dmult)` (the default, with
  `Dup=4`) exposes only the *certificate* degree: it takes PIETOOLS' `light`
  settings — whose Lyapunov degree is per-variable 1 — and sets the
  certificate degree to Lyapunov + `Dup` in every variable. So
  `set2d_deg(5,[])` enriches the certificate one degree further and, by
  construction, leaves the Lyapunov degree at 1. (The optional `dmult`
  raises only the Lyapunov block's multiplier cell — rarely needed.)
- **A full settings struct** changes anything, including the Lyapunov
  degree: e.g. `opts.settings = settings_PIETOOLS_heavy_2D()` (Lyapunov
  degree 2 per variable, with its own certificate degrees), or take
  `s = set2d_deg(4,[])` and edit `s.LF_deg` / `s.eq_deg` directly before
  passing it.

The certificate degree is why the eq Gram (`Ns(2)`, hundreds) dwarfs the LF
Gram (`Ns(1)`, single digits); it is usually the degree worth raising first
(§7, item 1).

## 7. When it does not certify

Work down this list — the failure signatures differ and mean different
things:

1. **The full problem is infeasible or beyond the LPI at these degrees.**
   Discovery residuals are O(1) at every rank, or the report says the
   reference solve failed. Raise the *certificate* degree first
   (`opts.settings = set2d_deg(5,[])`); if that is not enough, raise the
   *Lyapunov* degree by passing a full settings struct
   (`settings_PIETOOLS_heavy_2D()`) — see the note under §6. Or move the
   operating point away from the stability boundary: near the boundary of
   what the LPI can certify at all, minimum rank climbs toward full —
   borderline certificates are not compressible.
2. **Near-miss just above the gate** (residual ~1e-6–5e-6 at some rank):
   the search's characteristic false negative — a certificate at that rank
   may still exist. Add seeds, raise `.maxrank` by one or two, or supply a
   related face as a warm start. Accept the higher rank if it certifies;
   it is still a valid upper bound.
3. **Your operator forces a multiplier.** If the certificate must reproduce
   a multiplication operator `a(s1,s2)` exactly, its rank includes the
   number of squares `a` needs — unbounded in two variables — and a weight
   positive on the box but not a sum of squares is infeasible outright at
   these settings. See theory.pdf §7. Pure PDE stability problems dodge
   this structurally.
4. **Memory.** Assembly (not certification) is what grows: on the demo's
   heat-equation ladder — where `n` is the number of replicated PDE states,
   the `nlist` of `demo_lowrank_2d.m`, with Gram sizes `n × [8 744]` — it
   needs ~8 / 10 / 13 GiB at `n` = 4 / 5 / 6. Certification stays tiny at
   every size; the general rule is that memory follows the assembled
   program's size (constraint count grows like `n²`), not the face.

A failed run never produces a false certificate — that is the design — so
the cost of experimenting is only time.

## 8. Beyond stability: certifying a given operator

`pielr_certify_pos(Ptgt)` certifies **bare positivity** of a given `opvar2d`
operator — one Gram block, equated cell-by-cell to your `Ptgt` — with the same
option surface, chaining (`pielr_certify_pos(Ptgt, cert)`), acceptance gate
and report as the stability routine. `demo_poincare_2d.m` is the worked
example: the 2-D Poincaré inequality `∬a·u_{s1s2}² ≥ c·∬u²`, which unlike the
stability benchmarks is genuinely degree-hungry.

Two things change relative to stability:

- the certificate is **forced** to realise every cell of `Ptgt`, including a
  multiplier cell `R22{1,1} = a` if the target has one — its Gram sub-block
  then carries rank equal to the number of squares `a` needs, and a weight
  that is positive on the box but not a global sum of squares makes the
  problem infeasible outright (theory.pdf §7);
- the degree knob is `.spec = [n1 n2 n3]` (multiplier basis degree, then
  primary/dummy degrees of the integral bases; default `[2 2 1]`) — a
  starting spec, not derived from the target, so infeasibility at the default
  means raise it.

Tensor replication across states does not apply here; a previous certificate
still chains as a warm start on the same or a nearby target.
