# Ledger – D2 (PIESIM/, 62 files), 11/20/2026

Status: complete. Before the first edit, codeident.sh on all 62 files: 62 OK (no file skipped).
Final codeident.sh on all 62: **62 OK, 0 DIFF, 0 NEW**. Final runs.awk (N = 20) over all 62: **0 runs** (45 at
the start). Line endings: every file CRLF, final-newline state as at HEAD (re-derived from HEAD after the
scratchpad warning). 57 files differ from HEAD; fcgltran2d, fcgltran2dperm, fcht, ifcht, zwgll are identical to HEAD.
Whitespace: trailing blanks that an edit (11/19 or 11/20) removed from an otherwise unchanged line were restored to
HEAD, including 11/19 lines in PIESIM_1Dto1D2Mat_cheb_2D, PIESIM_1Dto2D2Mat_cheb_2D, PIESIM_2Dto1D2Mat_cheb_2D,
PIESIM_3PI_opmult_discretize, PIESIM_disturbance_check.

## Maintainer decisions of 11/20

- Support line: the 48 files in which the pass added `% For support, contact M. Peet, ...` now carry the two lines
  (`... mpeet@asu.edu,` / `% or Y. Peet at ypeet@asu.edu`): every D2 file except PIESIM.m, solver_PIESIM.m,
  rescalePIE, rescalePIE_2D, transl_1D, transl_2D, transl_3PI, transl_mod (HEAD lines kept) and the six third-party files.
- Third-party: the release log line removed from fcgltran2d, fcgltran2dperm, fcht, ifcht, zwgll by editing; each now
  matches HEAD (`git diff HEAD` empty). clencurt: no license, support or log line; keeps the 11/19 formula,
  INPUT/OUTPUT and step comments. fcgltran4d (YP) unaffected.

## Edited files (11/20)

Mechanical-only files (rules 1, 5, 6):
- PIESIM_Poly2Mat_cheb_2D: formula; N 1x2; Rop class; p corrected (2 x ns, header said scalar); output sizes; step comments.
- PIESIM_discretize_all_2D, PIESIM_discretize_all: PIE, uinput, psize and output structures; 1D psize fields corrected (no, n1, n2 are not fields).
- PIESIM_discretize_domain, PIESIM_discretize_domain_2D: grid formula (Chebyshev-Gauss-Lobatto nodes mapped to the domain); fields; sizes.
- PIESIM_discretize_ops, PIESIM_discretize_ops_2D: PIE, psize, Dop fields; 1D optional map fields named as set (Tchebmap_2PDEstate, ...); step comments.
- PIESIM_discretize_icf: formula; fields; coeff and B1_nonpol; MMP body comment restated without the pointer to the batch/legacy input checks (marker kept, column 77); step comments.
- PIESIM_fullPI2Mat_cheb_2D: block-matrix formula and truncation rule; Rop fields; psize; step comments.
- PIESIM_3PI2Mat_cheb: 3PI formula; Rop fields; output sizes; step comments.
- PIESIM_Poly2Mat_cheb: formula; one-variable restriction; p corrected (1 x ns, said scalar); A_2PDEstate size corrected ((n0+n1+n2)(N+1), said n0n1n2(N+1)^3); step comments.
- PIESIM_initial_setup: fundamental-state formula, 2D corner values, conversion of inputs to function handles (sym, constant, spline of samples); fields; step comments.
- PIESIM_input_check: header corrected (varargin = {structure, uinput, opts}, all required; header swapped uinput and opts, said 1-4 arguments); "checked in PIETOOLS converter" removed; psize/uinput/opts fields; "(pie2piesim)" removed from 3 MMP body comments (markers kept); INPUT/OUTPUT for both local functions; step comments (12 runs).
- PIESIM_options_check: defaults, resets and type detection stated; pointers to PIESIM_input_check and the converter removed; varargin -> varinput in the input list; help-block blank lines closed; step comments.
- PIESIM_int_bdf: BDF formula; input fields; solcoeff.final (header said acheb_f), timedep.coeff; release line moved to the true end of the log; step comments (3 runs).
- PIESIM_int_symbolic: convolution formula; fields; solcoeff.final; step comments.
- PIESIM_stability_check: boundary-locus criterion and thresholds; opts fields; step comments (62-line run).
- PIESIM_time_integrate: psize fields corrected (nf, no not fields); solcoeff.final; NOTE corrected (one other scheme, intScheme = 2, analytic; no intScheme 3 or Gauss scheme); timedep.coeff, w, u; step comments.
- examples_pde_library_PIESIM_1D, _2D: "(see below)" replaced by the case range (1-41, 1-32); output classes and uinput fields; help-block blank line closed.

Passed files (step comments and corrections):
- PIESIM.m: help-block blank lines closed; 2 time schemes (said 3, with a Gauss scheme); ifexact is an opts field (3 places); timedep.regulated indices {3}, {4} (said {2}, {3}); grid a 2x1 cell.
- PIESIM_9PI2Mat_cheb_2D: 9PI formula; Rop as 3x3 cell; prow/pcol 2 x n (said 1 x n); NOTE: 2PDEstate for T, Tu, Tw, Tmap (said "Mcheb"); step comments.
- PIESIM_9PI_opint_area_discretize_2D, PIESIM_discretize_icf_2D (formula, fields), PIESIM_3PI_opint_discretize, PIESIM_PI2Mat_opint_cheb: step comments.
- PIESIM_4PI2Mat_cheb: 4PI formula; Rop fields; step comments.
- rescalePIE, rescalePIE_2D: L, K and the maps listed as inputs; step comments.
- transl_1D: signature with var1, var2; R is one kernel; "default [-1,1]" removed (no default exists); step comments. transl_3PI: same default removed; step comments. transl_2D: T and I structures; step comments (50-line run).
- transl_mod: "(pie2piesim)" removed from a body comment (marker realigned to column 77).
- PIESIM_plot_solution, _2D: what is printed and plotted; psize fields corrected; uinput/grid/opts fields; 2D: observed size no0 (said n0), regulated indices; step comments (9 runs).
- PIESIM_transform_to_solution, _2D: formula; PIE input "not used" (header said needed); 1D psize fields corrected (nf); 2D Tcheb_nonsaquare -> Tcheb_2PDEstate, regulated indices; fields; step comments.

## Skipped

None.

## Candidates for the maintainer (read, not run; nothing fixed)

1. PIESIM_time_integrate: the BDF2 fallback (Atotal not diagonalizable) calls PIESIM_int_bdf with 5 arguments; it takes 6 (gridall missing), so Dop is undefined inside.
2. PIESIM_int_symbolic reads psize.no, which no file assigns (grep over the repository): intScheme = 2 stops there. It also indexes uinput.w(i), wdot(i), u(i), udot(i) as symbolic arrays, while PIESIM_initial_setup stores cells of function handles.
3. PIESIM_initial_setup (2D corner values): `if (length(uinput.w{i}==1))` and the u analogue take the length of a logical array, true for any nonempty input; the branch for 2 x m sampled inputs is never reached (`length(...)==1` intended?).
4. PIESIM_input_check: (a) PDE.dim == 0 calls initialize_PIETOOLS_PDE(psize,true) with psize an empty struct instead of the augmented PDE; (b) 2D PDE without uinput.ic uses undefined `n` (ns); (c) DDE branch: exact-solution check uses `ns`, undefined there; (d) PIE branch, 1D: control-input count tests ~isempty(PIE.w_tab) before reading PIE.u_tab; (e) legacy terms check compares the exact-solution count with psize.n0+psize.n (a vector).
5. PIESIM_transform_to_solution (1D): `if (ns>0)` in the NaN branch (ns undefined); `solution.timepde.pde=[]` (misspelled field). The PIE argument is unused in both transform files.
6. PIESIM_plot_solution, _2D: observed outputs 1-7 are labelled $z_i$ and colored by colors_y(psize.no0-ii+1,:), out of range for no0 > 7; 1D returns before assigning the output figs when opts.plot is not 'yes'.
7. PIESIM_9PI2Mat_cheb_2D: A_opmult and A_opint_block_* are assigned only for nonempty kernels but summed unconditionally (as 11/19 candidate 2 for 1Dto2D/2Dto1D).
8. PIESIM_discretize_icf (1D, w and u) and _2D (u): uinput.wspace(k)/uspace(k), 1x1 cells, passed to PIESIM_NonPoly2Mat_cheb, which applies subs to Rop(m,i); unverified that subs accepts the cell (2D w branch uses wspace{k}).
9. transl_mod translates the kernels but leaves T.I at the original interval.
10. Known, not fixed as instructed: PIESIM_discretize_icf_2D uses nu without defining it.

## Unresolved

1. PIESIM_9PI_opint_area_discretize_2D: sample order (x rows, theta columns, then y, nu) stated from the one-variable branches; the two-variable branch relies on subs of a row polynomial by a row vector (not run).

## For the maintainer

- clencurt keeps the 11/19 step comments besides the formula and INPUT/OUTPUT header (the decision names only the header).
- fcgltran4d (unchanged here): "Initial coding unknown" at HEAD, while its help says "Written by Yulia Peet 07/05/22".
- Literature citations kept: examples_pde_library_PIESIM_1D examples 35-36 ("Amritam's paper", other-author text, possessive kept), 38 (Lhachemi & Prieur 2021), 41 (A. Gahlawat et al., "Bootstrap Current Optimization in Tokamaks Using Sum-Of-Squares Polynomials"); clencurt (Trefethen, source URL).
- Change-log entries and (was) lines that name pie2piesim, pie2opvar or other files are kept (new = 0 files, rule 4).
