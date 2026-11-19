# A stability certificate 0.04% above the margin: PDE 57 (FO2s), D123 (MMP, 11/18/2026)

Report (Sister session, 08:30): in the B1 tightness bisections, candidate T2 (light, lpi.ineq.dD 0,
sop.deg_Q{2} 0 raised by degree twice to 2; 240,830 variables, m 7,821) certifies PDE 57 at
tau = 2.8825265 against the margin 2.881401, 0.04% above it, with MOSEK OPTIMAL, OPT_STATUS 0.867,
eta 1.5e-8, relres 8.6e-4. Log studies_1118/logs/T2_fo2s_1.log; record results_T2_fo2s.csv in
testing_codes/private/baselines/2026-11-18_studies/private. The maintainer is offline; this note
records what was measured and what is inferred, for his decision.

## 1. The system and its margin, independent of the solver
PDE 57 (PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman): z_tt = 0.1 z_ss - 2 z_t - z(t - tau, s)
on [0, pi], z = 0 at both ends. The sine modes are complete, so the PDE is stable iff every mode is:
mode k obeys y'' + 2y' + 0.1 k^2 y + y(t - tau) = 0 (the example file, MMP 11/17/2026). Mode 1:
lambda^2 + 2 lambda + 0.1 + e^(-lambda tau) = 0. Wolfram FindRoot at 20 digits: the crossing is at
omega = 0.4947333, tau* = 2.8814005 (library 2.881401); d Re(lambda)/d tau at the crossing 0.0608.
Rightmost real part of mode 1 at the bisection steps (the first value is the decay or growth rate):

| tau | Re lambda | harness verdict (T2) | MOSEK exit at the final rung | OPT_STATUS | relres | eta |
|---|---|---|---|---|---|---|
| 2.73733 | (stable) | C | OPTIMAL | 1.013 | 8.4e-4 | 1.0e-8 |
| 2.84538 | (stable) | C | OPTIMAL | 0.958 | 9.5e-4 | 6.2e-9 |
| 2.8724 | -5.50e-4 | C | OPTIMAL | 1.109 | 9.2e-4 | 5.3e-9 |
| 2.87915 | -1.37e-4 | C | OPTIMAL | 0.810 | 5.4e-4 | 9.4e-9 |
| 2.8814005 | 0 (tau*) | | | | | |
| 2.8825265 | +6.85e-5 | C (the reach) | OPTIMAL | 0.867 | 8.6e-4 | 1.5e-8 |
| 2.88421 | +1.70e-4 | U (harness gate) | OPTIMAL | 0.480 | 6.2e-4 | 1.9e-8 |
| 2.8859 | +2.73e-4 | F (harness gate) | ILL_POSED, PRIMAL_ILLPOSED_CER | -1 | 1.0 | 1.0 |
| 2.89941 | +1.09e-3 | F | UNKNOWN (divergent) | -1 | 7.7e7 | 1.8e-3 |
| 2.95344 | +4.25e-3 | F | UNKNOWN (divergent) | -1 | 6.7e6 | 1.7e-3 |
| 3.16954 | (unstable) | F | UNKNOWN (divergent) | -1 | 3.9e5 | 1.6e-3 |

MOSEK's own summary at every OPTIMAL exit of the table: primal constraint violation 2e-4 to 5e-4 at
a primal norm of 5e1 to 2e2 (relative 3e-6 to 6e-6), stable and unstable steps alike. No printed
quantity separates the two unstable OPTIMAL exits from the four stable ones: OPT_STATUS 0.867 and
0.480 against 0.810 to 1.109, relres 8.6e-4 and 6.2e-4 against 5.4e-4 to 9.5e-4, eta 1.5e-8 and
1.9e-8 against 5.3e-9 to 1.0e-8.

## 2. Two findings in the executive's own verdicts (read in the log; measured by the Sister's run)
(a) At tau = 2.88421 (growth 1.7e-4) the executive reports "optimal: certificate at epneg 0: 1".
Only the harness gate (|OPT_STATUS - 1| <= 0.3) refused it. The executive's acceptance reads the
solver flags alone (lpi_classify_sop: pinf, dinf, numerr; the eta gate of 11/05 was recommended
and not adopted), and eta is 1.9e-8 here, below every gate discussed.
(b) At tau = 2.8859 (growth 2.7e-4) MOSEK exits with problem status ILL_POSED and solution status
PRIMAL_ILLPOSED_CER, x = 0 (residual norm 1, eta 1). SOSSOLVE's MOSEK branch maps this exit to
pinf 0, dinf 0, numerr 0 (its tests read 'INFEASIBLE' and 'UNKNOWN' in the status strings;
'ILL_POSED' matches neither), lpi_classify_sop returns 'optimal', and the executive prints
"optimal: certificate at epneg 0: 1". The harness refused it through OPT_STATUS = -1 only. The
mapping lines (SOSTOOLS400/sossolve.m 404-421) are Declan's (2021, 2024): a defect report, not a
patch (Sec. 5).

## 3. Reading
- The certificate at 2.8825265 is a solver acceptance at a margin below its resolution, the family
  of the DEMO12 closed-loop bound and PDE 24 'Q' (memory classify-flags-only-negative-gap, 11/05); MEASURED in Sec. 4 (strict mode UNKNOWN at every tau):
  MOSEK's default INTPNT_CO_TOL_NEAR_REL 1000 reports OPTIMAL at a relative primal violation of
  about 5e-6, the same at every step of this program family, and the derivative constraint with
  epneg 0 claims no decay rate that a residual could be measured against (the maintainer's rule of
  09/28: convert the residual into an operator bound and compare it with the claimed margin; with
  epneg 0 the claimed margin is 0).
- On this system the solver's verdicts order by the growth rate: decay 1.4e-4 and growth 6.9e-5
  both OPTIMAL; growth 1.7e-4 OPTIMAL at OPT_STATUS 0.48; growth 2.7e-4 ILL_POSED; growth 1.1e-3
  and above infeasible by divergence. A reach bisection to 1e-3 of tau* on this program family is
  therefore not resolvable; the Sister reads reaches to 1e-2 (its rule), which is safe here.

## 4. Reproduction (ab_fo2s.m, Fable scratchpad libpass/ab; results_fo2s_margin.csv in the 11/17
## harness folder; log libpass/logs/ab_fo2s_rep.log)
Queued at 08:50 for the next MATLAB slot: tau in {2.8780, 2.8814, 2.8825265, 2.8842, 2.8859} with
the T2 program posed directly (light, dD 0, deg_Q{2} 2, no raise; the variable count checks the
identity with the Sister's final rung); each (1) through the executive with the default MOSEK
settings; (2) under MOSEK's strict near-optimal rule (INTPNT_CO_TOL_NEAR_REL 1) through a scratch
shadow of sossolve that also maps ILL_POSED to numerr 2; (3) at the default solve, the residual of
the derivative equality D + N = 0 in operator coefficients relative to max|P| and max|D|.
RESULTS (08:48-08:53, 319 s wall; the program is the Sister's final rung: 240,830 variables, m 7,821 at every tau):

| tau | Re lambda | default MOSEK: executive status, certificate | OPT_STATUS | relres | eta | MOSEK primal nrm, violation | strict near_rel 1 | operator residual max|D+N|, max|P|, max|D| |
|---|---|---|---|---|---|---|---|---|
| 2.878 | -2.07e-4 | optimal, 1 | 0.918 | 3.7e-4 | 1.0e-8 | 9e1, 1e-4 | UNKNOWN, 0 | 9.5e-8, 0.164, 0.0535 (rel_D 1.8e-6) |
| 2.8814 | -3e-8 | optimal, 1 | 0.629 | 4.0e-4 | 1.6e-8 | 5e1, 2e-4 | UNKNOWN, 0 | 1.6e-7, 0.0729, 0.070 (2.2e-6) |
| 2.8825265 | +6.85e-5 | optimal, 1 | 0.552 | 4.8e-4 | 1.7e-8 | 5e1, 2e-4 | UNKNOWN, 0 | 1.7e-7, 0.0685, 0.0732 (2.3e-6) |
| 2.8842 | +1.70e-4 | ILL_POSED, PRIMAL_ILLPOSED_CER: optimal, 1 | -1 | 1.0 | 1.0 | (x = 0) | UNKNOWN, 0 | 4e-13, 1.4e-12, 2.8e-12 (P = 0) |
| 2.8859 | +2.73e-4 | UNKNOWN (divergent): unknown, 0 | -1 | 1.1e6 | 9e-3 | 6e5, 3e5 | UNKNOWN, 0 | |

Read from the table:
- The Sister's certificate reproduces: at 2.8825265 the executive certifies (OPT_STATUS 0.552 here
  against 0.867 in its run; the harness gate would refuse this one). The ill-posed exit moved from
  2.8859 (Sister) to 2.8842 (here) and the divergent one from 2.89941 to 2.8859: the verdicts within
  0.2% of tau* are not reproducible to the step, a second sign of the resolution floor.
- Under MOSEK's strict near-optimal rule (INTPNT_CO_TOL_NEAR_REL 1) every program of the table is
  UNKNOWN, the stable one at decay 2.1e-4 included, at the same point (same OPT_STATUS, relres, eta).
  MEASURED: this program family is solved to a relative primal violation of 1e-6 to 4e-6 (MOSEK's
  own summary) and no further; every OPTIMAL exit on it is a near-optimal acceptance.
- The derivative equality D + N = 0 holds to 1.8e-6 to 2.3e-6 of max|D| in operator coefficients,
  stable and unstable alike: the certificate at 2.8825265 is as accurate as the one at 2.878.
- The ill-posed exit at 2.8842 reproduces finding (b) of Sec. 2 on this run: x = 0, P = 0 to 1e-12,
  "optimal: certificate at epneg 0: 1".

## 5. Defect report for Declan (SOSTOOLS400/sossolve.m, MOSEK branch, lines 404-421 at 96c9961b)
An ILL_POSED problem status (solution status PRIMAL_ILLPOSED_CER or DUAL_ILLPOSED_CER) is mapped to
pinf = 0, dinf = 0, numerr = 0, the values of an optimal exit, because the branch tests only for
'INFEASIBLE' and 'UNKNOWN' in res.sol.itr.prosta and res.sol.itr.solsta. The returned x is 0
(residual norm 1 after the b normalization). A caller that reads the flags takes a certificate of
ill-posedness for an optimal solution. Measured: the SOS program of PDE 57 at tau = 2.8859 (T2
configuration, 240,830 variables, m 7,821; log T2_fo2s_1.log lines 2225-2269). Suggested mapping:
numerr = 2 when prosta contains 'ILL_POSED' (the SeDuMi convention for a failed solve), so that
sossolve's callers classify the exit as unknown.

## 6. For the maintainer (decisions, nothing implemented)
- The ILL_POSED mapping: Declan's lines in SOSSOLVE (Sec. 5, report). Independently of that fix,
  lpi_classify_sop (MMP, 10/30/2026) could treat an exit whose residual is the whole right-hand side
  (relres >= 0.5 or eta >= 0.5, x = 0) as 'unknown' whatever the flags: a guard, not a tolerance.
- The near-optimal acceptance: no printed quantity separates the certified unstable steps from the
  certified stable ones on this program family (Sec. 1 and 4), so no gate on OPT_STATUS, relres or
  eta repairs it; a stricter MOSEK rule rejects the stable steps too. What separates them is the
  claimed margin: with epneg 0 the certificate claims no decay rate, and a residual of 2e-6 of max|D|
  admits a growth rate of the order measured here (6.9e-5 certified, 1.7e-4 not), the Trap C of
  09/28 (memory lpi-feasibility-acceptance-tests). A reach bisection on such a family is reliable
  to about 1e-3 of the margin at best; the Sister reads reaches to 1e-2.
- An example-level consequence: the library's "stable iff tau < 2.881401" for PDE 57 cannot be
  checked to better than 0.1% by this LPI at this degree; the st block of the example certifies the
  default 0.9 tau* with a margin of 2.6e-2 in decay rate, far above the floor.
