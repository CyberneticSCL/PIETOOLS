# Commit plan, 11/18/2026 (prepared on the maintainer's word; nothing committed yet)

## Basic verification before the commit
1. ab_sep.m (01:33): the six tiers parse and validate; the 11/17 programs reproduce with the
   fallback; PDE 1 and PDE 33 light unchanged; the new knobs act. PASSED.
2. verify_demos.m (queued 01:50): the tiers load; DEMO3, DEMO7, DEMO12 (Dan Braghini's 11/18 version)
   run end to end; the stability and H-infinity gain executives on PDE 1 and PDE 20 at light and
   heavy. PENDING.
3. Library pass: light and heavy complete and compared three ways (D91-D94); veryheavy to PDE 55;
   extreme and stripped running (part B); the DDE and NDSDDF light rows at deg_Q [4 2 2] queued.
4. Peer sessions: low-rank bit-identical to HEAD (30 fingerprints); cuADMM suites consistent,
   only the three light 1-D 'Q' rows moved (D95-D96 era messages).

## Proposed scope (one commit, "Settings tiers on the report layout; direct form separated")
- settings_sopvar/: the six tier files, lpisettings_sop.m, lpisettings_stock_sop.m, README.md,
  settings_2d_stock.m (new).
- executives_sopvar/: PIETOOLS_stability_sop.m, PIETOOLS_H2_norm_sop.m, PIETOOLS_Hinf_estimator_sop.m,
  PIETOOLS_Hinf_gain_sop.m, utility_functions/exec_sop_settings.m, exec_tools_sop.m.
- lpi_programming_sopvar/: lpi_ineq_sop.m, poscopvar.m, misc/get_lift_degs.m, misc/poscopvar_direct.m,
  misc/poscopvar_lift.m and misc/posmult_cdopvar.m (renames already staged from orphaned_file/).
- PIETOOLS_demos/: DEMO3, DEMO5, DEMO12 (merged with Dan Braghini's version); cuadmm/ and
  lowrank_2d_stability/ (the peer sessions' integration of the new settings, 9 files).
- PIETOOLS_examples/: examples_PDE_library_PIETOOLS.m (the [PDE, st] output, cases 53-72), the 2-D
  interior-delay example, Examples_Library/Mixed_Domain_Examples/ (new; the example-settings agent
  is adding st blocks there now: commit the state at commit time, or wait for its report).
- PIETOOLS_documentation/for_developers.txt.
- testing_codes/private/: the eight tracked files of the cuADMM, low-rank and mixed-domain tests;
  baselines/2026-11-17_tiers/ (new: the harness copy, results_A/B.csv, the HEAD reference run, the
  A/B scripts and logs).
- cleaning/private/settings_1116/ (new: decisions, run log, studies, the report and proposal with
  pdfs, about 160 files): removed before release; include unless the pdfs should stay out.

## Not in scope
- Values still under study: the tier tweaks (deg_Q by the rule, max_raise_Q, raise_Q 'deg') and
  the example-specific settings (the agent's EXAMPLE_SETTINGS.md).
- The converter session's leaner Sedumi2Mosek edit (held for the maintainer's answer).

## Commit message draft
    Settings tiers on the 11/16 report layout; direct stability form separated

    settings_sopvar: the six tiers carry a storage row per dimension (lengths,
    patterns, cap, cross length), lpi.ineq (dD, dw, wmin, wcross, psatz
    {'auto','faces1','auto'}), lpi.var.deg, sop (slack, dw_feas, dw_P, dw_Q,
    max_raise_feas, max_raise, cap_delta, margin, and since 11/18 deg_Q,
    max_raise_Q, raise_Q for the direct form 'Q'); settings_2d leaves the
    tiers for settings_2d_stock(level). Light's 1-D deg_Q is [4 2 2].
    executives_sopvar: the direct form 'Q' (the storage T'Q with Q free, two
    slacks sized from their own support, base raise 0) posed in every
    dimension; exec_sop_settings validates the new fields; exec_tools_sop
    reads the row in 2-D unless an explicit settings_2d is given.
    lpi_programming_sopvar: get_lift_degs gains wcross, dw_space and the
    term set 'faces1'; lpi_ineq_sop reads the tier's term set.
    Examples: optional second output [PDE, st] with application-specific
    settings; mixed-domain examples 56-72.
    Measured: library pass of 11/17 under all tiers against the 10/28
    reference and a HEAD reference run (cleaning/private/settings_1116,
    DECISIONS D84-D97); low-rank programs bit-identical; cuADMM suites
    consistent.

## Demo pass, 02:26-02:40 (verify_demos.sh: each demo its own MATLAB process)
- Tiers: all six load and validate (light deg_Q [4 2 2]).
- DEMO7 ok (35 s). DEMO3 failed at lpisolve(prob): sossolve refuses the empty solver field of
  the 11/17 tier layout -> lpisolve resolves '' like exec_sop_settings (staged); rerun pending.
- DEMO12 (Dan Braghini's 11/18 version): "for dw = 0:settings.sop.max_raise" with the per-dimension
  vector max_raise [1 1 0] -> "Colon operands must be real scalars" (his inlined loop; the
  executives reduce the vector in exec_tools_sop). Reported in the run log, not patched; the
  same on origin dd465ee5; the merged file is committed as is.
- The executive check of the first pass called the interactive examples library (not batch
  safe); rerun through the harness (ab_exec.m: PDE 1, 20, 20 dual, 33 at light and heavy).
