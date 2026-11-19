# DEMO12_Hinf_OFcontrol.m: failure report (MMP, 11/18/2026, for Dan Braghini)

Version: commits 78ab4501 and dd465ee5 (DB, 11/17/2026 19:16 and 20:38), merged into the
settings revision bb367fe1 on 11/18 with the 11/17 storage-row lines of the closed-loop
analysis kept. Run in batch mode on the toolbox at bb367fe1 (MATLAB R2025b, MOSEK).

## What fails
The PDE converts, then line 138

    for dw = 0:settings.sop.max_raise   % DB 11/17/2026

stops with "Colon operands must be real scalars". `settings` is `lpisettings_sop('light')`;
its `sop.max_raise` is a vector with one entry per number of spatial variables, [1 1 0] in
every tier since 10/29/2026 (settings_sopvar/README.md, "raises"; `max_raise_feas` likewise).
The executives reduce such vectors to the entry for the PIE's dimension in
`exec_tools_sop` (field list `dw_feas, dw_P, dw_Q, max_raise, max_raise_feas, cap_delta,
max_raise_Q`) after `exec_sop_settings` has filled the defaults; the inlined synthesis of the
demo reads the raw tier field. The same code fails on origin at dd465ee5, before the merge.

## Smallest fix (not applied here: your code)
Either reduce the fields the loop reads to the PIE's dimension, as the executives do,

    nv = PIE.dim;   mr = settings.sop.max_raise;   mr = mr(min(nv,end));
    for dw = 0:mr

or take the settings through the executives' two helpers once,

    st = exec_sop_settings(settings);   tl = exec_tools_sop(PIE,st);   st.sop = tl.sop;

and read `st.sop.max_raise`, `st.sop.dw_Q`, `st.sop.margin` (the last two are scalars in the
tiers today, but the same reduction protects them if a tier ever makes them per dimension).

## Two things to know about the settings this demo now depends on (since 11/17 and 11/18)
- `poslpivar_settings_sop(prog,T,settings,'lf','in',dom)` reads the storage row `lpi.pos`
  of the tier (lengths, patterns, cap); the former fields `mult, int, joint, sep` are refused.
- `lpi_ineq_sop(prog,-D,o)` with `o.like = R` sizes the slack like R; `o.dw` is the weight
  raise on top of the rule; the term set and the lift options come from `settings.lpi.ineq`
  of the program when the call gives none.
- The direct stability form reads its free operator from `sop.deg_Q` (11/18), not from
  `lpi.var.deg`; a hand-written Q form that copies the executive should do the same.

## Not checked
The PIESIM closed-loop part after the synthesis (never reached), and whether the inlined LPI
reproduces `PIETOOLS_Hinf_OFcontrol_sop` (the call is commented out at lines 121-122).

## Fix applied (maintainer, 11/18 03:35: "Fix the Demo, keeping the inline format")
Before the raise loop, `max_raise` and `dw_Q` are reduced to the entry for `PIE.dim` (the
first of the two smallest fixes above); the loop and its exit test read the reduced scalars,
and the slack options use the reduced `dw_Q`. Old lines kept as (was); header entry MMP
11/18/2026. The inlined synthesis and the PIESIM part are otherwise untouched. Run result:
see the run log entry of the demo pass (demo12.log).

## Second defect of the same kind (03:30)
With the loop fixed, line 159 failed: `poslpivar_settings_sop(prog,T,settings,'lf','in',dom)`
reads `settings.lpi.pos.int` (the stock-format storage fields mult, int, joint, sep), which the
tiers replaced by the storage row (len_inc, len_com, len_sos, ...) on 11/17/2026. The executives
build their storage through POS_ROW of EXEC_TOOLS_SOP: `poscopvar(prog,dm,sp,dom)` with DEG
omitted reads the row. Applied in the demo for R and V (`poscopvar(prog,dmi,spi,dom)`), the
old lines kept as (was). Anything else in the toolbox that still calls
POSLPIVAR_SETTINGS_SOP with tier settings fails the same way; EXEC_TOOLS_SOP keeps it only as
the fallback when no row is present.

## Run after the two fixes (03:32, 41 s)
Synthesis: gam = 5.6252653, certified at dw 0 (MOSEK OPTIMAL). Closed-loop analysis with
`PIETOOLS_Hinf_gain_sop(PIEcl,settings_an,'Q')`: 5.6236756, status optimal (residual 8.9e-5,
eta 1.5e-5), but the Legendre witness reports a numerical closed-loop gain of 5.6486 at
omega = 7.8e5, above both bounds, so the demo's guard of 11/05 prints "certifies no bound".
Either the witness responds to a high-frequency artefact of the truncated model or the two
certificates are loose by 0.4%: worth a look (N_cheb, nfreq of the witness; the analysis eta
is above the 1e-7 gate the library harness uses). PIESIM: the closed-loop integration completes
("numerically stable"), regulated outputs -0.0590 and 0.6359 at t = 1.998.
