# The container pipeline above 3 variables – 11/24/2026

Maintainer's question after the storage-row fix (NOTE.md): what else in the code does not
support more than 3 variables. Code survey of `lpi_programming_sopvar`, `executives_sopvar`,
`settings_sopvar`, `sopvar`, `converters` and `PIESIM`, then three MATLAB probes (R2025b,
HEAD 336b5fdc plus the uncommitted `poscopvar.m` fix) on the heat equation
x_t = sum_k x_{s_k s_k} + x, Dirichlet at both ends of every direction, in 3 and 4 variables.
Scripts and diaries in this folder: `probe4d_pipe.m` (`pipe_diary.log`), `probe3d_exec.m`
(`exec3d_diary.log`), `probe4d_strip.m` (`strip4d_diary.log`), runner `run_probe.sh` (memory
watchdog: the working set of the MATLAB processes of the run, polled every 10 s).

## Code: no cap at 3 anywhere else

Every per-dimension setting is read at its last entry for a higher dimension: `exec_tools_sop`
(dw_feas, dw_P, dw_Q, max_raise, max_raise_feas, cap_delta, max_raise_Q at `min(nv,end)`;
`settings_nd` at `min(3,end)` of lpi.pos int/mult/joint and lpi.var.deg), `lpi_ineq_sop`
(psatz cell, dD, dw, wmin at `min(max(numel(P.vars),1),end)`), `lpivar_cdopvar` (lpi.var.deg
at `min(max(a,1),end)`), `PIETOOLS_stability_sop` (deg_Q, dw_feas, max_raise_feas). The row
reader of `poscopvar` was the one exception (fixed, NOTE.md). `lpiprogram` accepts any number
of variables (the 2-variable refusal recorded 09/21 is gone). `eye_container` takes a scalar
margin or 2^nv entries. `get_lift_degs`, `poscopvar_direct`, `build_exponent_grid`,
`degbalance_core` and the cap of `lpi_ineq_sop` enumerate subsets as 2^n or 2^(2n) arrays and
kernel types as 3^n cells: generic in n, exponential in size.

Capped at 2 by design, not 3: PIESIM and `pie2piesim`; the GUI; the 'stock' slack and
`settings_2d` (`settings2possopvar`, `poslpivar_settings_2d_sop`); `getObserver_sop` (the
H-infinity estimator returns P and Z without a gain in 3 or more variables); the witness
(1-D). Capped at 1: the controllers, the H2 estimator, well-posedness.

## Measured

| Step | 3 variables | 4 variables |
|---|---|---|
| `pde_var` declaration + `convert_PIETOOLS_PDE` | 2.0 s | 2.5 s |
| `lpivar_cdopvar(T, deg_Q entry 3 = 1)` | 1,000 dvars, 0.1 s | 10,000 dvars, 0.0 s |
| light storage (`pos_row`) + margin | 4,225 dvars, 0.2 s | 66,049 dvars, 0.6 s |
| stripped storage | – | 4 dvars, 0.3 s |
| `T'*T`, `T'*A` (copvar) | – | 2.6 s, 2.4 s |
| D = T'PA + A'PT, light storage | 49–51 s, peak about 13 GB | > 31 GB at 157 s, STOPPED |
| `T'*P`, `(T'*P)*A`, stripped storage | – | 5.4 s; 133.5 s, peak about 18 GB |
| `lpi_ineq_sop` estimate, light (dw 1, cap Jmax(T) + 0 = 2 per direction) | 7 Gram blocks of 1,000: 7,000,000 certificate variables | 9 Gram blocks of 4,096: 150,994,944 |
| `lpi_ineq_sop` build | 85.8 s (light storage) | STOPPED at 26 GB, 384 s (stripped storage) |
| `PIETOOLS_stability_sop(PIE, light, 'P')` | SDP 7,004,225 variables, m 105,850; solve passed 33 GB at 287 s, STOPPED | not attempted |
| `PIESIM(PDE)` | error `Unrecognized field name "nx"` (PIESIM_input_check line 232) | – |

Stops are the watchdog of `run_probe.sh` (24 GB), not MATLAB errors. Machine: 64 GB.

## Reading

1. **Composition is the first wall in 4 variables.** `(T'*P)*A` with 4 decision variables
   takes 134 s and about 18 GB; with the light storage (66,049 decision variables) it passes
   31 GB and was stopped. `int_semisep` holds 3^n kernel cells per (beta, alpha) pair, 81 in
   4 variables against 27 in 3, each a polynomial in 8 variables. Measured, not profiled:
   which product dominates inside the composition was not localized.
2. **The 'like' slack is a tensor product per direction.** 10 monomials per direction at
   light give 1,000 per term in 3 variables and 4,096 (8 per direction under the cap) per
   term in 4, over 2n + 1 terms: 151M certificate variables in 4 variables. Nothing on record
   solves above 13.9M (Sol, 69 GB, 5 h). No tier value reaches a buildable 4-D slack;
   the lever is the form and the cap rule (slack-cap-sizing-nd memory), not a tier entry.
3. **3 variables at HEAD, light, 'P' form, workstation: not solvable.** 7.0M variables and
   m 105,850 on the Dirichlet heat equation; the solve was stopped at 33 GB. The 11/09 record
   (the N-D storage rule of `settings_nd`, before the storage row of 11/18) was 907,729
   variables and 289 s on the 3-D heat plant of that date; the 11/17 Sol record of the direct
   'Q' form at light is 4.07M variables, m 41,164, 39.7 GB. The storage row raises the
   storage Gram from the earlier rule to 65 (a plus b with 4 monomials per direction), and
   the 'like' slack follows the storage. The 3-D entries of the tiers were listed open in the
   11/16 settings report; this measurement does not close them. Different plants (Dirichlet
   here; mixed conditions in the library examples 53–55): the sizes are indicative, not a
   regression proof.
4. **PIESIM on a 3-D PDE fails on an internal field**, not on the dimension message of
   `solver_PIESIM` ("Dimension must be specified as 1 or 2"), because `PIESIM_input_check`
   reads `psize.nx` first. `PIETOOLS_PDE.m` calls PIESIM for every dimension (no gate before
   Step 4) and builds initial conditions for at most 2 variables (`nd = min(dim,2)`), so a
   3-D library example run through the entry script errors after its LPI. Inferred from the
   code of the entry script; PIESIM itself measured. PIESIM is YP's code: report, not patch.
5. **Tests.** `test_executives_sop` covers 1-D and 2-D; 3-D only in the baseline and study
   harnesses (heatNd); above 3 variables only `test_posmopvar` (09/21, class level) and these
   probes.

## Not done

No 4-D solve, no 4-D build completed, no profile of the composition, no 5-D declaration
(storage alone was declared in 5 variables, NOTE.md). Nothing edited in the toolbox for this
note.
