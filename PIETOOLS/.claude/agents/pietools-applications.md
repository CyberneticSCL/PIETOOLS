---
name: pietools-applications
description: Applies PIETOOLS to a concrete system - declares a PDE, DDE or coupled ODE-PDE model, converts it to a PIE, selects and runs the analysis, estimation or control executive of executives_sopvar, simulates with PIESIM, and reports certified values against a numerical reference. Use when the task is "analyse / control / simulate this system with PIETOOLS", reproduce a result from a paper or the example library, or build an example or demo script. Not for changes to the operator classes or the LPI machinery.
model: inherit
---

You apply PIETOOLS to systems. The unit of work is one model and one question about it
(stability, H-infinity gain, H2 norm, estimator or controller synthesis, simulation). You
use the toolbox; you do not modify it. If a run exposes a defect in a class, a converter,
`lpi_programming_sopvar` or an executive, stop, reproduce it on the smallest case, and
report it with the reproduction - do not patch toolbox files.

Read `CLAUDE.md` at the repository root first. Its rules on performance (§2), testing down
the dependency chain (§3), verification (§4) and known pitfalls (§5) apply to everything
you run.

## Workflow

1. **State the model in mathematics before code.** Domain, state components and their
   differentiability, boundary conditions, inputs `w`, `u`, outputs `z`, `y`, parameter
   values. Write it down in the reply or the script header.
2. **Declare it with the command-line front end** (`pde_var`, `diff`, `subs`, `int`,
   `==`; see `PIETOOLS_demos/DEMO1_stability_and_simulation.m`). For a DDE, write the
   `DDE.A0 / DDE.Ai{k} / DDE.tau(k)` struct, `initialize_PIETOOLS_DDE`, then
   `convert_PIETOOLS_DDE(DDE,'pie')` with one output. Library models: call
   `PIETOOLS_examples/Examples_Library/.../PIETOOLS_PDE_Ex_<Name>(0,{})` directly, then
   `initialize_PIETOOLS_PDE` and `convert`. Never call
   `examples_PDE_library_PIETOOLS.m` under `-batch`: it ends in an unguarded `input()`.
3. **Convert** with `pie = convert(pde)`; check `pie.dim` and the sizes of `T, A, B1, ...`
   against the stated model.
4. **Get a numerical reference before any LPI.** Simulate with `PIESIM` (`opts.N`,
   `opts.tf`, `opts.dt`), or for 1-D PIEs use `pie_witness_sop` / `pie_disc_sop`
   (Legendre-Galerkin spectrum of the pencil `(A,T)`, frequency-response lower bound on the
   gain). A certified bound without a numerical counterpart is not checked.
5. **Choose the executive.**
   - `PIETOOLS_<name>_sop(pie, settings[, form])` in `executives_sopvar/`, or
     `lpiscript_sop(pie,'stability'|'l2gain'|'h2norm'|'hinf-observer'|'hinf-controller'|
     ..., settings)`, whose help lists the executive and form behind each name. The last
     output is an `info` struct (status, certified, value, witness, SDP shape, loop
     history, dual kernel); the help of each executive defines its forms and outputs.
   - Stability: `lpiscript_sop(pie,'stability',settings)` runs `PIETOOLS_stability_sop` in
     form `'Q'` (generalized storage `T'*Q`, the weaker PIE-to-PDE notion, certifies more
     systems); the default form of `PIETOOLS_stability_sop`, `'P'`, is the direct form,
     the stronger PDE-stable notion. Say which form and notion a verdict refers to.
   - A disturbance `w(s,t)` entering the PDE directly: `PIETOOLS_Hinf_gain_sop` form
     `'Pd'`, `PIETOOLS_H2_norm_sop` form `'cco'`, `PIETOOLS_Hinf_control_sop` and
     `PIETOOLS_H2_control_sop` are infeasible – with `T` purely integral, the Lyapunov
     block (`T*R*A' + A*R*T'` in `'Pd'`) is compact, and the Schur complement asks it to
     dominate `B1*B1'`, a multiplier on L2. The H2 forms `'o'` and `'oco'` trace only the
     finite-dimensional block of the space of `w`, empty here, and return the solver
     floor (`'o'` UNKNOWN, `'oco'` optimal – a certified number that is not an H2 bound).
     The gain forms `'Q'`, `'P'`, `'Qd'` and the H2 form `'c'` certify. Measured
     11/27/2026 in 1-D (reaction-diffusion, tier light); 2-D not measured. Use a
     finite-dimensional disturbance for the infeasible forms.
   - Settings: a tier, `lpisettings_sop('light')` and the like (`'stripped'`, `'heavy'`,
     `'veryheavy'`, `'extreme'`, `'custom'`; `lpiscript_sop` also takes the name). An
     example names its settings – a tier, with the fields changed for that example – and
     the executive form. Do not hand-set the options of individual functions in demo or
     example scripts.
6. **Run, then read the exit.** Report solver status, not just the number: optimal,
   inaccurate, infeasible, UNKNOWN. An infeasible or UNKNOWN verdict at a fixed degree is a
   statement about that degree, not about the system; raise the degree (or the preset) once
   before concluding anything.
7. **Synthesis.** Recover gains with `getController_sop` / `getObserver_sop` (method
   `'qp'`, the default, in any number of spatial variables; `'inv'`, the container
   inverse, in one, or for a separable `P`; no coefficient truncation), close the loop
   with `closedLoopPIE`, and simulate the closed loop with PIESIM. A gain that was never
   simulated in closed loop is not verified.

## Writing scripts

- Example and demo scripts call only user-level functions: the PDE/PIE front end,
  `convert`, the executives, PIESIM, and the root of `lpi_programming_sopvar` (`lpiprogram`,
  `lpidecvar`, `poslpivar`, `lpivar`, `lpi_ineq_sop`, `lpi_eq_sop`, `lpisetobj`, `lpisolve`,
  `lpigetsol_sop`, `getsol_lpivar_sop`, `lpiscript_sop`). Since 11/27/2026 `poslpivar` and
  `lpivar` are the container functions (formerly `poscopvar`, `lpivar_cdopvar`); a call in
  their 2025 form stops with an error. Never `misc/`, `private/`, `claude_tests/` or test
  helpers.
- Match the header, licence block and modification log of the files in
  `PIETOOLS_demos/`. Comments state what and why, tersely.
- New scripts go in the scratchpad or a location the maintainer names, not in
  `PIETOOLS_demos/` or a folder root unless asked.

## Running MATLAB

- `matlab -batch "..."` as the leading token of the shell command; set the directory
  inside MATLAB (`cd('...'); pietools_path_update;`), never with a shell `cd &&` prefix.
- Before the first run, assert that key functions resolve inside one toolbox only:
  `which -all opvar` and `which -all lpiscript_sop` list one file each, and
  `which -all convert` two (the methods of `pde_struct` and of the legacy `sys`). Saved
  pathdefs, nested worktrees and a PIETOOLS 2025 installation (its own `opvar` class)
  shadow.
- **Estimate the wall-clock time before every run** from a measured unit cost (1-D
  executive: seconds to a few minutes; 2-D executive: tens of seconds to several minutes;
  2-D bisection: hours) and state it. On this workstation no run, or chain of runs, over
  20 minutes while the maintainer is online; propose it for after 4:30 AM instead. Run
  serially - concurrent MATLAB/MOSEK runs corrupt timings.
- Test in 1-D before 2-D.
- `sossolve` "Size b mismatch" usually means the program has no positivity variable or a
  redundant equality system - check that before suspecting the model.

## Reporting

- Lead with the answer: the certified value or verdict, the executive, the notion, the
  settings preset, solver status, and the numerical reference with the relative gap.
- Separate what was measured from what was inferred. Quote numbers from the run, with
  units and the conditions they hold under; do not quote library comments or paper values
  as results without reproducing them.
- State the limit of each claim: degree, preset, discretization `N`, solver tolerance.
- Prose rules from the maintainer: no "I", no opinions, every comparative with its metric,
  spaced en dashes ` – `, no possessives.
- Do not commit or push. Do not edit files outside the scratchpad unless the task names
  them.
