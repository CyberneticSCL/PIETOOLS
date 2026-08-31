# PIETOOLS — instructions for Claude Code

Conventions agreed with the maintainer. These are requirements, not suggestions: the
licence block of every source file already mandates the documentation rules below.

---

## 1. Code modification protocol

Applies to every edit of an existing file.

**Header entry.** Add to the file header: developer initials, date, and a brief summary of
the change *and the reasoning for it*. Follow the format the file already uses — these vary
and you should match the file you are in rather than impose one style:

```matlab
% Initial coding MMP, SS, DJ  - 09/26/2021      % lpi_programming/poslpivar.m
% DJ, 10/16/2024: Update to new 'lpiprogram' structure;
% 05/04/23 - DJ: Allow slightly negative eigenvalues.   % SOSTOOLS400/custom/findsos.m
```

**Per-line markers.** Append your initials and date to each changed line, corresponding to
the header note. Align the `%` at **column 77**, matching the existing markers throughout
the toolbox:

```matlab
    degs_i = unique(degmat_full(:,idcs1_new(i)));                           % XY, 08/24/2026
```

Do not include line numbers — they go stale as soon as anything above shifts. Never place a
marker on a line ending in `...`; it belongs after the continuation.

**Minor changes — keep the old code, commented out**, stamped with date and initials:

```matlab
%   maxdeg = max(Pjj.degmat(:,var));                                        % XY, 08/24/2026 (was)
    maxdeg = max(Pjj.degmat(:,vidx));                                       % XY, 08/24/2026
```

**Extensive changes — mark the start and end** of the modified section with initials and
date. Where commenting out the old code is impractical, state in a comment what was deleted
and the date stamps of the code being removed. If the change invalidates an earlier header
entry, note that at the end of the header log so a reader knows the earlier change is no
longer present or no longer relevant.

**Say what and why, everywhere.** Every portion of the body — minor change, extensive
change, or new code — states what is being done and why, citing the governing mathematical
documentation for the class or the user manual where appropriate. This applies especially to
modifications and exceptionally to major ones.

**Comment wording: terse.** Coverage is non-negotiable, wording is not. Explain everything,
in as few words as carry the fact. Cut `Note that…` lead-ins, restatements of what the code
plainly does, narration of how a bug was found, and rationale for roads not taken. Keep the
non-obvious invariant, the reason a wrong-looking line is correct, and the failure mode a
guard prevents. Do not restyle comments written by other authors.

---

## 2. Performance is a correctness requirement

**Typical use involves millions of decision variables across several spatial dimensions.**
Evaluate every proposed change for time *and* memory at that scale before offering it.

In `ndopvar`/`nopvar`, `q = numel(Pop.dvarname)` is the millions-scale quantity, and it enters
the coefficient matrices as the **row** dimension:

```
Pop.C{i}  is  dim(1)*nZ*(q+1)  x  dim(2)*nZ_t
```

Rows of `C` are therefore the dangerous axis. The monomial and spatial axes (`deg`, `nZ`,
number of directions, quadrature nodes) are small; dense work there is fine and common.

- Never `full()`, `zeros()`, `ones()` or `repmat()` on a dimension involving `q`. Use
  `sparse(m,n)` — it allocates only the column-pointer array, so an empty sparse with
  millions of rows is nearly free.
- Build sparse matrices once from triplets: accumulate `i`, `j`, `v`, then a single
  `sparse(i,j,v,m,n)`. Never grow a matrix by concatenation or index-assignment inside a
  loop over decision variables. Note that a transient `[A; B]` of `q`-row blocks doubles
  peak memory even when the result is sparse.
- Prefer remapping nonzero indices over `full` → `reshape` → `permute` → `reshape`.
- Skip empty or all-zero blocks rather than filling them. Fill only where dimensions
  genuinely cannot be inferred otherwise.
- `unique`, `sort`, `setdiff` and `find` on `q`-length arrays are each a hidden
  `O(q log q)`. On the spatial-variable axis (a handful of entries) they are harmless.
- **When a hot helper conflicts with a rare caller, fix the caller.** Do not slow a routine
  that runs on every decision variable to accommodate an edge case in a debugging
  convenience.

**Cost is multiplicative across spatial directions.** `nZ = prod(deg+1)` and the number of
parameter cells both grow with the number of spatial variables, so a change that is flat in
`q` at one spatial dimension can still blow up at three or four. Sweep a grid over (number
of spatial variables, `q`) rather than testing one large-`q` point.

**Measure, do not reason.** A static argument that a change is `O(1)` in `q` is a hypothesis.
Time it and profile memory over a sweep. State each fix's cost in `q` and in spatial
dimension, for time and memory, alongside the correctness argument. If the only correct fix
is `q`-scaled, say so explicitly rather than quietly shipping it.

---

## 3. Testing scope and dependencies

Test as extensively as possible: **varying numbers of spatial variables**, and replicating
the intended workflow of solving large-scale SDPs for analysis, control and simulation of
non-trivial PDEs.

```
pvar / dpvar  ->  opvar / dopvar / nopvar / ndopvar / sopvar / sdopvar
opvar / dopvar  ->  lpi_programming 1D + 2D  ->  executives  ->  examples
nopvar / ndopvar  ->  polyopvar / lpis_ndopvar
polyopvar  ->  tensopvar / intop  ->  sosprogramming
```

**Changes to core data structures must be tested all the way down this chain.** A green
class-level unit test on `dpvar` or `opvar` is not evidence the change is safe; exercise it
through `lpi_programming` and an actual executive or example.

Class locations (one definition each; `objs/opvar/` is inert):

| Class | Path |
|---|---|
| `pvar` / `polynomial` | `SOSTOOLS400/multipoly/@polynomial` |
| `dpvar` | `SOSTOOLS400/dpvar/@dpvar` |
| `opvar`, `dopvar` | `opvar/@opvar`, `opvar/@dopvar` |
| `opvar2d`, `dopvar2d` | `opvar/2D/@opvar2d`, `opvar/2D/@dopvar2d` |
| `nopvar`, `ndopvar` | `ndopvar/@nopvar`, `ndopvar/@ndopvar` |
| `sopvar`, `sdopvar`, `mopvar` | `sopvar/@sopvar`, `sopvar/@sdopvar`, `sopvar/@mopvar` |
| `polyopvar`, `tensopvar`, `tensopmat`, `intop` | `polyopvar/@…` |

`pietools_path_update.m` is a single `addpath(genpath(...))`, so everything is on the path at
once — there is no module isolation to rely on.

---

## 4. Verification discipline

**Never test mutually inverse routines against each other.** `ndopvar2sdopvar` and
`sdopvar2ndopvar` once mapped the decision coefficient with the same index error in exactly
inverse ways: the round trip returned the original object bit for bit while the intermediate
object was a *different operator*. A round-trip test is structurally incapable of detecting a
pair of inverse errors. The same pattern produced the `sopvar` ZL/ZR swap. Test each
direction against the semantics — evaluate the kernels from each class's own definition.

**Reproduce a reported defect before acting on it.** If a report hedges ("may", "could be"),
treat it as a hypothesis to test, not a finding. Distinguish in writing between what was
measured and what was inferred.

---

## 5. Known pitfalls

Each of these was expensive to discover. Dates are given because the code moves — verify
before relying on a specific line.

**The `dpvar` classdef header is wrong.** For `dpvar(C,degmat,varname,dvarname,matdim)`, rows
are *matrix row outer, decision variable inner* — one block of `nd+1` consecutive rows per
matrix row, constant term first, equivalently a left factor `kron(I_m,[1 d'])`. Columns are
matrix column outer, monomial inner. The classdef header states the decision variable is
outer; that is incorrect (verified 2026-08-28).

**`sopvar` ZL/ZR.** `ZL` holds monomial degrees in the *output* variables `P.vars.out`, `ZR`
those in `P.vars.in`. The two conventions agree whenever `numel(ZL)==numel(ZR)`, which is the
common case — with equal sizes but different degree *sets* a swap returns a silently wrong
answer with no error. Guard by checking `numel(ZL)` against `numel(vars.out)` and the
parameter shape. Any test using equal-size bases cannot catch a swap.

**Canonical multiplier form is a class invariant** (since 2026-08-29). Where the multi-index
has `gamma_k == 1`, `delta(s_k - s_k')` identifies the variables, so stored coefficients are
not determined by the operator; both constructors enforce a canonical form. *The adjoint is
the subtlety*: in an integral direction `ZL` and `ZR` swap; in a multiplier direction they
must **not**, because delta collapses the kernel to the diagonal. Monomial bases must be
**column** vectors — MATLAB's `union` returns a row when both arguments are rows, so always
append `(:)` to a `union` result used as a basis.

**`sossolve` "Size b mismatch" is usually not your bug.** It fails inside SeDuMi's
`pretransfo` for a program with only equality constraints and no SOS or psd variable
declared, and for some highly redundant equality systems — both reproducible with stock
PIETOOLS. Check whether the program has a positivity variable at all, and confirm whether the
stock routine fails identically, before suspecting new operator code.

---

## 6. Reference documentation

The mathematical spec for the `sopvar` / `sdopvar` classes is maintained on Overleaf
(*ndopvar Structure / alt_nopvar*, `sopvar.pdf`) and is **not** in this repository — ask the
maintainer for access. Section map: §1–3 cell-of-parameters representation, addition,
adjoint; §4–5 the kernel form and its adjoint; §6 composition; §7 `dsopvar` and the four key
operations; §8 the decision-variable object `vec C(d) = A + B'd`; §9 construction of positive
operators. Cite it by section when explaining a modification (protocol item 5).
