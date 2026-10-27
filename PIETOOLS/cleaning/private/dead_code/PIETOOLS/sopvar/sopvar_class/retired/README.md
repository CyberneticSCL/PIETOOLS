# Retired code

Superseded implementations, kept for reference. **Nothing here is on the
MATLAB path**: the folder sits inside an `@`-folder, and `genpath` excludes
class folders and everything beneath them. Confirmed by inspecting
`genpath(PIETOOLS)` — this directory does not appear. Files here can therefore
never shadow a live function, and are reached only by opening them directly.

This mirrors the existing `@quadPoly/incompleteFunctions` convention.

| file | was | removed |
|---|---|---|
| `int_semisep_local_mtimes.m` | local subfunction in `@sopvar/mtimes.m`, lines 666-1015 | 2026-08-30 |
| `rearrangeCoef_old.m` | local subfunction in `@sopvar/mtimes.m`, lines 1016-1105 | 2026-08-30 |
| `int_semisep_AT_v1.m` | `@sopvar/private/int_semisep_AT.m` | 2026-08-30 |

## Why they went

`int_semisep` existed in four near-identical copies:

- `sopvar/int_semisep.m` — called by `possopvar`
- a local subfunction in `@sopvar/mtimes.m`, which **shadowed** the shared file
- a byte-identical local subfunction in `@sopvar/mtimes_AT.m`, which was
  **dead** (`mtimes_AT` calls `int_semisep_AT`, not `int_semisep`)
- `@sopvar/private/int_semisep_AT.m`

A fix or an optimization applied to one silently missed the others. The lookup
table that replaced `ismember(...,'rows')` in the key loop had gone into the
shared file alone, so `@sopvar/mtimes` — the main composition path — never
received it.

They are now one implementation: `sopvar/int_semisep.m`, called directly by
`possopvar` and `@sopvar/mtimes`, and through the repacking wrapper
`@sopvar/private/int_semisep_AT.m` by `@sopvar/mtimes_AT`.

## The one behavioural difference

The copies disagreed on **repeated rows** in `idxbeta`/`idxalpha`: the shared
routine filled only the first occurrence, the local copies filled every row.
Both callers pass `unique(...,'rows','stable')`, so no live path could reach
it, but the shared routine now fills every occurrence so that the result does
not depend on callers deduplicating.

Equivalence was established before removal: 144 exhaustive index combinations
and random multi-variable cases at exactly 0 difference, both implementations
matched against the defining integral identity by quadrature, and `mtimes_AT`
driven end to end through the new wrapper agreeing with `mtimes` to 1.8e-15.

`rearrangeCoef_old` had no callers anywhere in the repo before any of this.
