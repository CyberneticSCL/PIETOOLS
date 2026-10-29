# lpi_programming_sopvar — LPI programming for both operator families

Initial coding MMP, 09/29/2026 (Tier 1 of the container parity map).
MMP, 10/06/2026: the 1-D executive translators (Tier 2); moved 10/16/2026, see below.
MMP, 10/16/2026: reorganized. This folder holds the LPI programming functions only. The
executives' helpers moved to `executives_sopvar/utility_functions/`, the constructors to
`sopvar/constructors/`, `op2copvar_sop` to `converters/operator_converters/`, `on_registry_sop` to
`sopvar/misc/`, and the slack-sizing studies to `sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/`
(§3). The measurements and the history that were in this file are in `NOTES.md`, verbatim.
MMP, 10/16/2026: `sopvar/lpis_sopvar` merged into this folder: the container routines the
dispatchers call (second table below) and their 11 `private/` helpers. Its settings translator
`settings2possopvar` went to `executives_sopvar/utility_functions/`.
MMP, 10/21/2026: the folder root holds the user functions only, the seven container counterparts
of the user functions of `lpi_programming/` (first table). Every routine they depend on is in
`misc/` (second table), and the 11 helpers of those routines in `misc/private/`.
MMP, 10/21/2026: `examples/` holds four scripts that use the user functions only (§4). The four
end-to-end checks that were there, which call test helpers, moved to `tests/` as `e2e_*`.
MMP, 10/21/2026: `lpi_ineq_sop` counts and prints the decision variables of its certificate before
declaring it, and sizes it from the operators of `poscopvar` its target depends on (§2); the
examples no longer set `opts.like`.
MMP, 10/28/2026: `lpiscript_sop`, the container counterpart of `lpiscript`, which had none
(maintainer: create one); the folder root holds eight user functions.

This folder sits beside `lpi_programming/` and lets one LPI body serve both operator
families:

- the legacy classes: `opvar`, `dopvar`, `opvar2d`, `dopvar2d`, `dpvar`, `polynomial`
  and `double`;
- the container family: `sopvar`, `sdopvar`, `copvar` and `cdopvar`, in any number of
  spatial variables.

Every `_sop` function sends a legacy class to the existing legacy routine, unchanged,
so its result is the legacy result. It sends a container-family class to container
code. No file in `lpi_programming/` was modified, apart from `lpiprogram`, which takes
any number of variables since 10/01/2026.

The user functions, in the folder root. Each is the container counterpart of the user function of
`lpi_programming/` in the second column; `lpidecvar`, `lpisetobj` and `lpisolve` are used unchanged.

| function | stock counterpart | what it does |
|---|---|---|
| `lpiprogram_sop` | `lpiprogram` | declares the program; also from a container's registry (§2) |
| `lpivar_cdopvar` | `lpivar` | declares a free `cdopvar` decision operator between concatenations of mixed spaces |
| `poscopvar` | `poslpivar` | declares a positive self-adjoint `cdopvar` on a concatenation of mixed spaces; calls `poscopvar_direct` (since 10/17/2026; `copquadvar` before, the same operator and decision variables); records it in `prog.sop_pos` for `lpi_ineq_sop` (10/21/2026) |
| `lpi_eq_sop` | `lpi_eq` | imposes `P == 0` (§2) |
| `lpi_ineq_sop` | `lpi_ineq` | imposes `P >= 0`: sizes the positive slack, prints its decision variables, declares it, then `P - Pop == 0` (§2) |
| `lpigetsol_sop` | `lpigetsol` | the solved value of a decision variable or operator (§2) |
| `getsol_lpivar_sop` | `getsol_lpivar` | the same for a decision operator (§2) |
| `lpiscript_sop` | `lpiscript` | solves a predefined LPI by name (`'stability'`, `'l2gain'`, `'h2norm'`, `'hinf-controller'`, ...) with the container executive and form of `PIETOOLS_auto_execute_sop`; the outputs of `lpiscript`, then the executive's `info` |

The routines the user functions depend on, in `misc/` (public, since the executives, the tests and
the demos call some of them). From `sopvar/lpis_sopvar` (10/16/2026) and the former folder root
(10/21/2026); each can still be called directly:

| function | what it does |
|---|---|
| `subs_dvar_sop` | `lpigetsol_sop` at given decision values, for certifying a candidate point (§2) |
| `lpigetdual_sop` | the solver's multipliers of one container equality as its signed dual kernel; needs `prog.sopeq = {}` set after `lpiprogram` |
| `lpi_shape_sop` | the SDP shape: decision variables, PSD block sizes, equality rows, primal length |
| `lpi_classify_sop` | the solver exit: optimal, inaccurate, infeasible, unbounded or unknown |
| `copquadvar` | the quadratic form in a list of basis operators, with one Gram matrix over every (space, basis operator) pair, built from the pair integrals; types `'pos'` and `'sym'`; the reference of the tests of `poscopvar_direct` |
| `poscopvar_direct` | the operator of `poscopvar`, mapped from the Gram blocks straight to the kernel coefficients; used by `lpi_ineq_sop` |
| `poscopvar_lift`, `lift_copvar`, `posmult_cdopvar` | the separated form `Pop = A'*M*A`: the fixed moment lift `A`, the positive multiplier `M` |
| `possopvar`, `sopquadvar` | the one-space (`sdopvar`) counterparts of `poscopvar` and `copquadvar`; `sopquadvar` calls `copquadvar` |
| `get_lift_degs` | sizes the positive operator that cancels a self-adjoint `P`, for `lpi_ineq_sop` |
| `eq_opts_sopvar` | the `possopvar` options read off the support of an operator (counterpart of `get_eq_opts_2D`) |
| `lpi_eq_cdopvar`, `lpi_eq_sdopvar` | impose `P == 0` for a `cdopvar` / `sdopvar`; called by `lpi_eq_sop` |
| `lpi_soseq` | appends to `prog.expr` the entry `soseq` would write for affine scalar equations, built from triplets |

## 1. Naming rule

`pietools_path_update` adds every folder to the path with one `genpath`, so a file here
must never reuse a legacy name. Such a file would shadow the legacy one for every
caller. Every dispatcher here therefore ends in **`_sop`**: `lpiprogram_sop`,
`lpi_eq_sop`, `lpigetsol_sop`, and so on. The container routines name their class
(`lpivar_cdopvar`, `poscopvar`, `possopvar`, `lpi_eq_cdopvar`, ...) or use a name no legacy
file has (`lpi_soseq`, `get_lift_degs`, `lift_copvar`). Helpers follow the same rule:

- `tests/private/same_val_sop` and `tests/private/opcheck_sop` are reachable only
  from `tests/` (`opcheck_sop` was in `examples/private/` until 10/21/2026, §4);
- `misc/private/` holds the 11 helpers of the container routines (`collect_eq_rows`,
  `impose_eq_rows`, `build_exponent_grid`, `unpack_sheets`, ...), reachable only from
  `misc/`, where every routine that calls one is. No other file in the tree has one of their
  names, so neither the 10/16/2026 merge nor the 10/21/2026 move changed the resolution of a
  call. A routine moved to `misc/` stays on the path (`genpath` adds `misc/`, never a
  `private/` folder), so every caller outside this folder still finds it.

MMP, 09/30/2026: `private/spaces2meta_sop` (a verbatim copy of the space parser of
`lpivar_cdopvar`) is deleted. The constructors call the shared
`parse_copvar_spaces` in `sopvar/misc/conventions/` instead, a name that does not
collide with a legacy one.

`which -all` resolves each function to one file in the tree.

## 2. Dispatch tables

### `prog = lpiprogram_sop(...)`

This is `lpiprogram` with three more input forms (the last three rows below). Since
10/01/2026 `lpiprogram` itself takes any number of variables, so `lpiprogram_sop` only
converts those forms and calls it; before, it was a copy of lpiprogram without the
2-variable cap. Its outputs and error messages were checked unchanged on 93 input forms.

| input | result |
|---|---|
| every `lpiprogram` form: `(vartab,dumvartab,dom[,decvartab[,freevartab]])`, `(vartab,dom,...)`, an n×2 `vartab`, a 1×2 `dom` for all variables | as `lpiprogram`, for any n |
| `(names,dom,...)` with `names` a cellstr or string, e.g. a container's registry `P.vars`, `P.dom` | polynomial variables in that order; dummies `s_dum` |
| `(P[,decvartab[,freevartab]])` with `P` a `copvar`, `cdopvar`, `sopvar` or `sdopvar` | `P`'s sorted registry and domains; dummies `s_dum` |
| `(struct('vars',..,'dom',..),...)` | the same |

It has these fields: `vartable = [vars; dummies; free]` (a `polynomial`), `dom` (n×2),
`decvartable` (a cellstr), and the SOSTOOLS fields.

The former copy deviated from `lpiprogram` in three places. Since 10/01/2026 `lpiprogram`
has two of them itself (it checks the dummy variables with `ispvar`, and with no spatial
variable sets `dom = zeros(0,2)`); the cellstr and container forms stay here.

### `prog = lpi_eq_sop(prog,P[,opts])`

| class of `P` | routine |
|---|---|
| `cdopvar`, `copvar` | `lpi_eq_cdopvar` |
| `sdopvar`, `sopvar` | `lpi_eq_sdopvar`, which collects and imposes the rows (MMP, 09/30/2026: it has one output and no internal modes since then) |
| any other class (`dpvar`, `dopvar`, `dopvar2d`, `opvar`, `opvar2d`, `polynomial`, `double`, ...) | `lpi_eq`, unchanged |

`opts` (`'symmetric'`) is passed on only when given. Each routine keeps its own errors,
for example for a fixed `copvar`/`sopvar`, a polynomial, or an unknown option. The
program it returns equals the routine's own program field by field (test_lpi_eq_sop).

### `[prog,Pop,info] = lpi_ineq_sop(prog,P[,opts])` (MMP, 10/08/2026)

| class of `P` | routine |
|---|---|
| `cdopvar`, `copvar` (square, self-adjoint; a `copvar` is certified) | this file: `get_lift_degs` sizes the positive operator from `P`, ONE `poscopvar_direct` call declares every Positivstellensatz term, `lpi_eq_sop` imposes `P - Pop == 0` ('symmetric') |
| any other class | `lpi_ineq`, unchanged |

The degrees come from the lower-kernel support of `P` per direction (`Dmin`, `Dmax`, the
multiplier degree): lift `D = Dmin + dD`, weight `w = max(ceil((Dmax - D - 1)/2), ceil(Mdeg/2), 0) +
dw`, defaults `dD = dw = 1` (the heat benchmark needs both, `PIETOOLS_demos/sopvar_demos/README.md`
Sec. 12); terms `'auto'`: 1-D the plain term at `w` and the product at `w - 1` (the Markov-Lukacs
pair), N-D the plain term at `w` and the `2N` faces at `w - 1` when every weight is at least 2, at
`w` otherwise (changed 10/08/2026, see below); `opts.deg` bypasses the reader; `opts.like` reads
the support of another operator instead, or of a cell of them, pooled. Without `opts.like` and
`opts.deg` it reads the operators of `poscopvar` that `P` depends on, which `poscopvar` records in
`prog.sop_pos`: a Lyapunov or KYP inequality is sized from its positive decision operator, not from
every value of it (10/21/2026, §4). Before anything is declared, the decision variables of the
certificate are counted, `N^2` per Gram block as `prog.decvartable` counts them, and printed with the
program total (`info.ndec`, `info.N`); `opts.prune`
(default true) keeps only the basis operators that can reach the diagonal blocks' support
(`eq_opts_sopvar`). `tests/test_lpi_ineq_sop`: a fixed positive operator certified, the 1-D heat LPI
feasible at kappa = 9 and rejected at 11, 2-D `T'*T` certified, the option variants, the legacy
dispatch. The several-term call of `poscopvar_direct` is itself checked against the sum of
single-term calls, through `poscopvar_direct` and through `poscopvar` + `plus_batch`
(`claude_tests/test_poscopvar_direct` part 6: same decision variables, coefficients to rounding).

The measurements behind these defaults (the H-infinity gain LPI in 1-D and 2-D, the coercive and
dual forms, the Volterra norm, the term sets) are in `NOTES.md`.

### `sol = lpigetsol_sop(prog,X)`, `Psol = getsol_lpivar_sop(prog,P)`, `Xsol = subs_dvar_sop(X,names,vals)`

These were written by the getsol agent. They are verified again here in the chain and
through E1–E4.

| class | `lpigetsol_sop` | `getsol_lpivar_sop` |
|---|---|---|
| `sdopvar` → `sopvar`, `cdopvar` → `copvar` | new code | new code |
| `sopvar`, `copvar` | returned as is | returned as is |
| a cell holding any of the four classes | taken apart, one element at a time | — |
| `[]`, `double`, `polynomial`, `opvar`, `opvar2d`, `dpvar`, char, another cell, `dopvar`, `dopvar2d` | `lpigetsol`, unchanged | `opvar`…`dopvar2d`: `getsol_lpivar` |
| anything else | error `lpigetsol_sop:badClass` | error `getsol_lpivar_sop:badClass` |

Values are read from `prog.solinfo.RRx`, never from `solinfo.x`: RRx is in decvartable
order, and x is the solver's cone vector.

- Each operator's decision list `Zd` is looked up **once** per container.
- The solved parameter is `A_γ + B_γ'd` for each γ cell. No `dpvar` of q names is
  formed.
- A name that is listed twice takes its value from its lowest row, the row where
  getequation places it.
- `subs_dvar_sop` is the same code applied to a given (names, values) pair, for
  certifying a candidate point.

### dpvar as an operator (Tier 1b, in the class files)

These changes are in the class files, not in this folder: `sopvar/@*/mtimes`, `plus`,
`horzcat` and `vertcat`, plus `sopvar/misc/dpvar_op_copvar.m`, `mat2copvar_grid.m` and
`dpvar2sdvar.m`. Each class method has one branch at the top of its body, and the
composition code below it is untouched.

They follow the legacy `opvar`/`dopvar` semantics:

- `gam*X` and `X*gam` for a scalar dpvar;
- `D*X` for a dpvar matrix;
- `gam ± X`, read as `gam*I ± X`;
- a dpvar entry inside `[ ]`, e.g. `[-gam, D'; D, -gam]`.

A product of two decision operators is refused (`cdopvar:decisionTimesDecision`).
Numeric operands keep their earlier behaviour: `P + 1e-3` and `[X, 0]` are still
refused. Use `dpvar(c)` or `mat2copvar_sop` instead.

### Legacy routines that need no `_sop` version

`lpidecvar`, `lpisetobj`, `lpisolve` and the scalar `lpi_ineq(prog,gam)` (sosineq) work
unchanged on container programs. **Measured** in two places:

- E1: 1-D, gamma declared, `gam >= 0`, and minimized;
- test_lpiprogram_sop (d): N = 1, 2, 3, 4 programs from `lpiprogram_sop`, with
  `lpidecvar`, `lpi_ineq`, `lpisetobj`, `lpisolve` and `lpigetsol_sop`.

Neither `prog.vartable` nor `prog.dom` is read by lpidecvar, lpisetobj, lpisolve or any
container routine.

## 3. Related functions in other folders (moved 10/16/2026)

- **Constructors**, `sopvar/constructors/`: `mat2copvar_sop`, `eye_copvar_sop`,
  `zeros_copvar_sop` (below).
- **Conversion**, `converters/operator_converters/` (in the PIETOOLS root): `op2copvar_sop` (an
  `opvar` or `opvar2d` to a container) and every converter between the operator classes that was
  in `sopvar/Testfolder/converters/`: `opvar2copvar`, `opvar2d2sopvar`, `opvar2d2copvar`,
  `opvar2sopvar`, `mat2sopvar`, `nopvar2copvar`, `nopvar2sopvar`, `ndopvar2sdopvar`,
  `dopvar2sdopvar`, `sdopvar2dopvar`, `dpvar2sdvar`, `dpvars2sdvars`, `sdvar2dpvar`.
- **Registry**, `sopvar/misc/on_registry_sop`: restates a container on a larger variable registry,
  which `@cdopvar/plus` needs when the operands' registries differ.
- **The container executives' helpers**, `executives_sopvar/utility_functions/`: the settings and
  degree loop, the stock settings translators (`poslpivar_sop`, `poslpivar_settings_sop`,
  `poslpivar_settings_2d_sop`, `settings_2d_sop`, and `settings2possopvar` from
  `sopvar/lpis_sopvar`), `get_lpivar_degs_sop`, `trace_rn_sop`, the gain
  reconstruction (`getController_sop`, `getObserver_sop`, `getController_direct_sop`) and the
  numerical witness (`pie_disc_sop`, `pie_witness_sop`); see `executives_sopvar/README.md`.
- **Slack-sizing studies**, `sopvar/Testfolder/sdopvar/claude_tests/slack_sizing/`:
  `hinf_tailor_1d`, `hinf_tailor_2d`, `volterra_tailor` (results in `NOTES.md`).

### Constructors (Tier 1c, in `sopvar/constructors/`)

The arguments `dims`, `spaces` and `dom` follow the conventions of `lpivar_cdopvar` and
`poscopvar`:

- `dims` is M×1, or `struct('out',..,'in',..)`;
- `spaces` is a 1×M cell of cellstr, where `{}` is R^q, or `struct('out',..,'in',..)`;
- `dom` is nv×2 in the order of the sorted registry, one 1×2 row for all variables, or
  `struct('vars',..,'dom',..)`.

| function | result |
|---|---|
| `mat2copvar_sop(Mat,dims,spaces,dom[,opts])` | the container of a constant matrix: `copvar` for a `double`, `cdopvar` for a `dpvar`. A block acts as a multiplier in shared variables, is constant in output-only variables and is integrated over input-only ones. `opts.mult_only` refuses an integral block, as `mat2opvar` does. |
| `eye_copvar_sop(dims,spaces,dom)` | the identity. Square only (`eye_copvar_sop:notSquare`). |
| `zeros_copvar_sop(dims,spaces,dom)` | the zero operator. Rectangular is allowed. |

All three pass `verify` and satisfy the canonical multiplier form. They replace the
`opvar2copvar(mat2opvar(...))` route, which has no N-D form.

## 4. Running

Set the path with `run pietools_path_update`, then remove the `.claude`/`worktrees`
entries. A saved pathdef or a nested worktree can shadow the tree. The solvers
are MOSEK 11 and SeDuMi (E3).

### Examples (MMP, 10/21/2026)

The scripts in `examples/` show standard usage, as the DEMOs of `PIETOOLS_demos/` do. A PDE is
declared with `pde_var` and converted with `convert` (or an `opvar` is declared), its operators
become containers with `op2copvar_sop`, and one LPI is written with the user functions of the folder
root and the stock `lpidecvar`, `lpisetobj` and `lpisolve`. Each runs by name, with MOSEK when
`mosekopt` is on the path and SeDuMi otherwise. Measured with MOSEK 11 on the development
workstation, the time including `convert`:

| script | LPI | result | decision variables | time |
|---|---|---|---|---|
| `volterra_norm_sop` | DEMO2: min `gam`, `gam*I - T'*T >= 0` | 0.636804 (2/pi = 0.636620) | 126 | 1 s |
| `stability_1d_sop` | `P >= 1e-3*I`, `T'*P*A + A'*P*T <= 0`, `x_t = x_ss + (pi^2/2)x` | feasible | 2041 | 2 s |
| `hinf_gain_1d_sop` | the KYP LPI of `PIETOOLS_Hinf_gain`, then `PIETOOLS_Hinf_gain_sop` | 0.033311, both | 2274 | 2 s |
| `stability_2d_sop` | the stability LPI of `x_t = x_s1s1 + x_s2s2 + pi^2 x` on `[0,1]^2` | feasible | 197920 | 19 s |

In the three LPIs linear in the decision operator `P`, `lpi_ineq_sop` sizes the certificate from `P`,
which `poscopvar` recorded in the program (§2), and prints its decision variables before declaring
it; no example sets an option. Sized from the target instead, the certificate covers every `P` of the
declared degree: 9225 and
27426 decision variables in 1-D, with the same verdict and gain, and 12036880 in 2-D, unsolved after
900 s and 31 GB. With `opts.dw = 2`, the executives' default for the coercive forms, the 2-D LPI has
761680 decision variables and takes 50 s; `PIETOOLS_stability_sop` (form 'P') certifies the same PDE with
its own storage in 6 s (187984 decision variables).

No example has three spatial variables: `convert` reduces a PDE to at most two. The 3-D LPI is
checked by E4 below.

### Tests

The end-to-end checks E1–E4, `tests/e2e_*` (the examples of this folder until 10/21/2026), compare
the container LPIs with the legacy ones. They use the cx_exec helpers
(`sopvar/Testfolder/sdopvar/claude_tests/cx_exec`), `heatNd_apply` (`PIETOOLS_demos/sopvar_demos`),
both on the path from `genpath`, and `tests/private/opcheck_sop`.

| command | what | time |
|---|---|---|
| `test_lpiprogram_sop` | (a) 16 lpiprogram forms field by field and 8 error messages; (b) N = 3, 4 against the hand-built N-D program; (c) the 7 new input forms; (d) N-D solves (Volterra N = 1, 2; multiplier N = 1..4) | 92 s |
| `test_lpi_eq_sop` | dispatch of 16 cases field by field, the option pass-through, 7 errors | 2 s |
| `test_getsol_sop`, `bench_getsol_sop` | getsol (getsol agent) | 18 s, minutes |
| `test_dpvar_ops_sop`, `test_dpvar_hinf_chain_sop` | dpvar operators and constructors (operators agent) | 17 s, 10 s |
| `test_endtoend_sop` | E1 `e2e_hinf_gain_1d` and E2 `e2e_stability_1d` by default; `'all'` adds E3 `e2e_stability_2d` and E4 `e2e_heat3d_stability` | 15 s; E3 5 min; E4 15–20 min |
| `test_translators_sop` | the executive translators: legacy dispatch, errors, the settings pair, Q degrees against stock, trace against its definition, space lists against opvar dims (MMP, 10/06/2026) | 9 s |

## 5. Measurements and history

`NOTES.md`: the slack sizing of `lpi_ineq_sop` measured on the H-infinity and Volterra LPIs, the
Tier 1 scope as of 09/29/2026, the costs, the end-to-end results E1–E4 and the 09/29/2026
regression chain, each as first written here.
