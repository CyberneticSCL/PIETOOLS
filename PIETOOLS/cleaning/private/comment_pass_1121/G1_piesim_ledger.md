# G1_piesim ledger – comment pass of 11/21/2026 (42 files under PIESIM/)

Per file: class, TAG, the note verbatim, apply counts, internal notes removed, open questions.
Checks per file: cpass check PASS, codeident OK, `git ls-files --eol` unchanged, numstat = numstat -w.

## Summary
- Classes: condense 42 (all edited); asnew 0; review 0; deferred 0. Skipped: none (no file had edits
  of another session at its first edit).
- TAGs: `YP, MMP` 33, `MMP` 5 (comment-only changes), `YP` 2 (third-party fcgltran2d, fcgltran2dperm),
  `YP, DJ, MMP` 1 (PIESIM.m), `VJ, YP, MMP` 1 (PIESIM_transform_to_solution.m).
- apply totals: dropped 16 (was) lines (development steps), kept 17 (was) lines (2025 code),
  retagged 46 markers; 4 BEGIN/END lines retagged by hand (9PI2Mat_cheb_2D, PIESIM.m).
- Final sweep (fixed cpass.pl): check PASS 41, file 18 two false positives (reworded 2025 entries with
  underscore dates); codeident OK 42; EOL and numstat unchanged 42; the modified files under PIESIM/
  are exactly these 42. Parse check (syntax_check via mlrun): 42 files, 0 with parse problems.
- Open questions for the maintainer:
  1. TAG authorship: DJ kept in PIESIM.m and VJ in PIESIM_transform_to_solution.m (real post-2025
     changes, since superseded); DJ left out of PIESIM_input_check.m (only a revert-of-merge commit);
     DJ/VJ left out where their initials occur only in reworded 2025 entries (10, 40, 41); YP kept in
     transl_2D.m although YP added and removed the same line (net nil).
  2. A global rename no -> n0 after 2025 corrupted words on 16 lines of 5 files (n0w, n0t, den0tes,
     polyn0mial, en0ugh, n0te), comments and 2025 entries, and one user message (code, not fixed):
     'Warning: solution final time does n0t match user input final time' (PIESIM_transform_to_solution.m).
  3. Commented "data array (double)" / "symbolic" input alternatives in the example libraries kept as
     user-facing alternatives (not development steps); remove if not wanted.
  4. Stale or inaccurate text outside the scope, left: 2025 comments on rsize/csize (01); PIESIM.m help
     lists ic.pde/ic.pie/ic.ode; YP entry "(2) pass options to examples" of solver_PIESIM.m had no code.
  5. Read, not run: PIESIM_Poly2Mat_cheb_2D.m grows its padded coefficient array past (N_x+1)x(N_y+1)
     when deg+2 > N_x+1 or N_y+1.

## 01 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_9PI2Mat_cheb_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
  %   renamed A_2PDEstate, input p renamed prow; N is [N_x, N_y], and prow,
  %   pcol are 2 x (number of states) arrays [p_x; p_y]; each block is built
  %   on the full grid of prod(N+1) coefficients, then its columns are
  %   truncated by pcol (A_2PDEstate) and its rows by prow (A); a term with
  %   an empty kernel adds zero.
  ```
- replaced entries: MMP 11/19 (release pass), MMP 11/20 (terms initialized to zero); YP code
  changes had no entry. BEGIN/END pair retagged by hand.
- apply: dropped 0, kept 0 (was), retagged 0.
- internal notes removed: none.
- open: 2025 comment lines "rsize ... / csize ..." describe variables no longer in the code (left,
  2025 lines).

## 02 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_NonPoly2Mat_cheb_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: N is [N_x, N_y], each
  %   state block has prod(N-p+1) rows; the function line carries the file
  %   name (was PIESIM_NonPOLy2Mat_cheb_2D); help and comments revised.
  ```
- replaced entries: MMP 11/19 (release pass). apply: dropped 0, kept 0, retagged 0.
- internal notes removed: none. open: none.

## 03 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_Poly2Mat_cheb_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
  %   renamed A_2PDEstate (ns*prod(N+1) rows, not truncated); N is [N_x, N_y]
  %   and p is a 2 x ns array [p_x; p_y]; R is evaluated on deg+2
  %   Chebyshev-Gauss-Lobatto nodes per direction (deg the maximal degree of
  %   R), its coefficients padded to (N_x+1) x (N_y+1); A is A_2PDEstate with
  %   the rows of each state truncated by p.
  ```
- replaced entries: YP 1/25/2026 (non-square matrix for the 2-D state reconstruction), MMP 11/19.
- apply: dropped 0, kept 0, retagged 0.
- internal notes removed: the empty body marker line `% YP 1/25/2026` (no content).
- open: (1) the 2025 marker `% DJ, 12/16/2024` on `for j=1:length(acheb)` was deleted after 2025
  (the line keeps trailing blanks); not restored. (2) read, not run: when deg+2 > N_x+1 or N_y+1
  the padded array grows past (N_x+1) x (N_y+1) and Alocal past prod(N+1) rows.

## 04 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_discretize_all_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: third output grid
  %   replaced by gridall, the physical grids of the states of every degree of
  %   differentiability; help and comments revised.
  ```
- replaced entries: MMP 11/19. apply: dropped 0, kept 0, retagged 0. internal notes removed: none.

## 05 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_discretize_domain_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output grid removed;
  %   N is [N_x, N_y], and gridall.x{i}, gridall.y{i} hold N_x-i+2 and N_y-i+2
  %   nodes; help and comments revised.
  ```
- replaced entries: `% 2/12/2026: YP - modified grids handling ...` (2 lines; a format prep does not
  recognize as an entry), MMP 11/19. apply: dropped 0, kept 0, retagged 0.
- internal notes removed: none.

## 06 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_discretize_icf_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: initial conditions read
  %   from the single list uinput.ic (ODE states first, then the states of x,
  %   of y, of (x,y)), not from uinput.ic.PDE and uinput.ic.ODE; psize.n0
  %   replaces psize.no; N is [N_x, N_y] and each 2-D state has its own
  %   [p_x; p_y]; uinput.wspace is a cell, and a disturbance that does not
  %   separate (uinput.wsep{k} false) takes w at t=0 as its spatial part;
  %   control-input coefficients: loops over nuy, nu2 (were nwy, nw2), row
  %   count prod(N+1)*nu2 (was (N+1)^2^nu2), nu assigned (was unassigned).
  ```
- replaced entries: YP 1/7/2026, YP 2/17/2026, MMP 11/19, MMP 11/20.
- apply: dropped 0, kept 0, retagged 1 (nu=psize.nu).
- internal notes removed: none.

## 07 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_discretize_ops_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output fields renamed –
  %   Mcheb_inv to Tcheb_inv, Mcheb_nonsquare, Tucheb_nonsquare,
  %   Twcheb_nonsquare to Tcheb_2PDEstate, Tucheb_2PDEstate, Twcheb_2PDEstate;
  %   the optional map is PIE.misc.Tmap, output Tchebmap_2PDEstate (was
  %   PIE.misc.T0, Mcheb0_nonsquare); C1, C2 discretized with flag 1 and
  %   Tu, Tw, B1, B2 with flag 2 (were 2 and 1); help and comments revised.
  ```
- replaced entries: MMP 11/19 (the 2025 entry YP 06/26/2025 stays). apply: 0/0/0.
- internal notes removed: none.

## 08 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_fullPI2Mat_cheb_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
  %   renamed A_2PDEstate; flag 1 now marks state-to-output and flag 2
  %   input-to-state operators (swapped); outputs may be infinite-dimensional
  %   – every row block is discretized when its kernel is nonempty, output
  %   rows with degree of smoothness zero; N is [N_x, N_y] and each 2-D state
  %   has its own [p_x; p_y].
  ```
- replaced entries: YP 1/27/2026 (infinite-dimensional outputs), MMP 11/19. apply: 0/0/0.
- internal notes removed: none.

## 09 PIESIM/PIESIM_routines/2D/Discretization_Of_Operators/PIESIM_integral_projection1D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: input csize and output
  %   Abig removed; A is the square (N+1)x(N+1) matrix (was its first csize
  %   columns); help and comments revised.
  ```
- replaced entries: YP 1/31/2026, MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 10 PIESIM/PIESIM_routines/2D/Postprocessing/PIESIM_plot_solution_2D.m
- class condense; TAG `YP, MMP, 11/21/2026` (prep proposed DJ too: its line `DJ, 07/17/2024 ...` is the
  2025 entry with "." changed to ";", so it stays and DJ is not a post-2025 author)
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: solution fields ode, pde
  %   replaced by primary{1,2,3,4} (ODE states; states of s1, of s2, of
  %   (s1,s2)), and observed, regulated are cells with the finite-dimensional
  %   outputs in {1}; psize fields n0, nr0, no0 (were no, nro, noo); grid is a
  %   2x1 cell of the grids in s1 and s2; opts.ifexact replaces
  %   uinput.ifexact; at tf = 0 the final states are plotted (was no plot);
  %   titles show solution.tf; error plot with color bar, font size 22.
  ```
- replaced entries: YP 01/07/2026, YP 02/17/2026 (3 lines), MMP 11/19. apply: 0/0/0.
- internal notes removed: none.

## 11 PIESIM/PIESIM_routines/2D/Postprocessing/PIESIM_transform_to_solution_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: input grid renamed
  %   gridall; solution fields ode, pde replaced by primary{1,2,3,4} (ODE
  %   states; states of s1, of s2, of (s1,s2)), and observed, regulated are
  %   cells {1,2,3,4} likewise (outputs of s1, s2, (s1,s2) in physical space);
  %   uinput.w, uinput.u are cells, a non-separable input has its spatial
  %   coefficients recomputed at each time; Dop fields *_2PDEstate (or
  %   *map_2PDEstate) replace *_nonsquare; psize.n0 replaces psize.no; N is
  %   [N_x, N_y].
  ```
- replaced entries: YP 12/31/2025, YP 1/27/2026, YP 02/17/2026 (3 lines), MMP 11/19. apply: 0/0/0.
- internal notes removed: none.

## 12 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_3PI2Mat_cheb.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
  %   renamed A_2PDEstate; the multiplier and integral terms are discretized
  %   as (N+1)x(N+1) matrices, then truncated to the resolved columns
  %   (A_2PDEstate) and rows (A); p, pcol may be row or column vectors.
  ```
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 13 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_4PI2Mat_cheb.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
  %   renamed A_2PDEstate; help and comments revised.
  ```
- replaced entries: MMP 11/19 (the 2025 entry `YP - ... 6_2_2025`, "observe" -> "observed", stays).
  apply: 0/0/0. internal notes removed: none.

## 14 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_NonPoly2Mat_cheb.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: input gridall renamed
  %   grid; help and comments revised.
  ```
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 15 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_Poly2Mat_cheb.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output A_nonsquare
  %   renamed A_2PDEstate; p may be a row or column vector; help and comments
  %   revised.
  ```
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 16 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_discretize_all.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: third output grid
  %   replaced by gridall, the grids of the states of every degree of
  %   differentiability; help and comments revised.
  ```
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 17 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_discretize_domain.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output grid removed;
  %   help and comments revised.
  ```
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 18 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_discretize_icf.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: initial conditions read
  %   from the single list uinput.ic (ODE states first, then the PDE states),
  %   not from uinput.ic.PDE and uinput.ic.ODE; psize.n0 replaces psize.no;
  %   no non-polynomial forcing when nw = 0 (B1_nonpol = []).
  ```
- replaced entries: YP 1/7/2026, YP 2/17/2026, MMP 11/05/2026 (4 lines), MMP 11/19.
- apply: dropped 0, kept 1 (was) (2025 code `if isfield(uinput,'Bpw_nonpol')`), retagged 3.
- internal notes removed: none (the two marked comment lines on the nw > 0 guard state the failure
  mode it prevents; kept).
- check: 2 failures, both false positives – lines 57 and 59 are the 2025 entries `YP 6_16_2022 ...`
  and `YP - ... - 6_1_2025`, reformatted by YP on 02/18/2026 to `YP, 6/16/2022:` and `YP, 6/1/2025:`;
  the 2025 wording uses underscore dates, which give cpass no 2025 key. Kept as they are (brief: a
  reworded 2025 entry stays). codeident OK.

## 19 PIESIM/PIESIM_routines/Discretization_Of_Operators/PIESIM_discretize_ops.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: output fields renamed –
  %   Mcheb_inv to Tcheb_inv, Mcheb_nonsquare to Tcheb_2PDEstate; new fields
  %   Tucheb_2PDEstate, Twcheb_2PDEstate; the optional maps are PIE.misc.Tmap,
  %   Tumap, Twmap, outputs Tchebmap_2PDEstate, Tuchebmap_2PDEstate,
  %   Twchebmap_2PDEstate (were PIE.misc.T0, Tu, Tw, outputs
  %   Mcheb0_nonsquare, Tucheb0, Twcheb0); help and comments revised.
  ```
- replaced entries: YP 1/6/2026 (2 lines), MMP 11/19 (the 2025 entry `YP - 6/16/2022`, continuation
  reworded T0 -> Tmap by YP, stays). apply: 0/0/0. internal notes removed: none.

## 20 PIESIM/PIESIM_routines/PIESIM.m
- class condense; TAG `YP, DJ, MMP, 11/21/2026` (DJ: pie2opvar conversion of 10/29/2026, since
  replaced by MMP; DJ is a post-2025 author of the file)
- note:
  ```
  % YP, DJ, MMP, 11/21/2026: Changes since PIETOOLS 2025: output time removed;
  %   grid is the primary-state grid, an (N+1)x1 array in 1D, a 2x1 cell
  %   {x grid; y grid} in 2D; a 'pie_struct' PIE is converted to structs of
  %   its kernels (no operator class); infinite-dimensional inputs and outputs
  %   are no longer an error; a field uinput.weight = k multiplies T, Tu, Tw
  %   by s^k (new local function weightPIE; error with ODE states; maps in
  %   PIE.misc.Tmap, Tumap, Twmap; was uinput.weighted, with T, Tw, B2 and
  %   PIE.misc.T0, Tu, Tw); in 2D the weight is ignored with a warning.
  ```
- replaced entries: YP 1/22/2026, DJ 10/29/2026 (2 lines), MMP 11/04/2026 (4), MMP 11/19,
  MMP 11/20/2026 (4, incl. "maintainer, 11/20/2026" quote).
- apply: dropped 3 (was) (the `pie2opvar` lines, DJ 10/29/2026 code, not 2025), kept 0, retagged 6.
- BEGIN/END of the local function weightPIE: 2-line BEGIN ("moved unchanged from
  PIESIM_Setup/weightPIE.m (one caller ...)") reduced to `% TAG: BEGIN local function weightPIE`;
  END retagged.
- internal notes removed: the history in the BEGIN line (above). Kept: `% Initial coding: YP 1/7/2026`
  in the help of weightPIE (Initial coding line).
- open: help text lists uinput.ic sub-fields ic.pde/ic.pie/ic.ode, while the setup routines read
  one list uinput.ic (help text not in scope; left).

## 21 PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_initial_setup.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: new input PIE and output
  %   psize; initial conditions read from the single list uinput.ic (ODE
  %   states first), not uinput.ic.PDE; the temporal part of a separable input
  %   becomes a function handle of t (sym through matlabFunction, constant,
  %   or 2 x m samples [t; w] through a cubic spline, zero outside the sample
  %   times), a non-separable input stays sym; in 2D, corner values of
  %   boundary inputs become finite-dimensional inputs placed first (psize.nw,
  %   nw0, nu, nu0 updated) and the spatial parts of the inputs are
  %   differentiated as the states they enter.
  ```
- replaced entries: YP 12/31/2025 (2 lines), YP 1/22/2026 (2), YP 2/6/2026, YP 2/17/2026, MMP 11/19,
  MMP 11/20/2026 (2).
- apply: dropped 4 (was) (YP code of 2026 replaced by the MMP 11/20 fix: `length(w{i}==1)`, scaling of
  row 1), kept 0, retagged 4.
- internal notes removed: none.

## 22 PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_input_check.m
- class condense; TAG `YP, MMP, 11/21/2026` (prep proposed DJ too: the only DJ commit since 2025 is
  48905ac4, a revert of a merge that changed one disp text, rewritten since; no DJ content remains)
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: new output opts
  %   (opts.ifexact replaces uinput.ifexact; opts.N scalar in 1D, 1x2 in 2D);
  %   psize fields n0, nr0, no0 replace no, nro, noo; new fields nw0, nwx,
  %   nwy, nw2 (likewise nu, nr, no) count the inputs and outputs by domain;
  %   in 2D psize.n is the matrix by orders in x and y; initial conditions in
  %   the single list uinput.ic, missing ones zero (default zero, was
  %   sinusoidal), extra ones dropped; inputs are checked and split into
  %   temporal and spatial parts; DDE and PIE sizes are read from field dim.
  ```
- replaced entries: YP 12/29/2025, YP 1/27/2026, YP 2/17/2026, MMP 11/04/2026 (two, 8 lines),
  MMP 11/19, MMP 11/20/2026 (2).
- apply: dropped 9 (was) (YP code of 2026), kept 6 (was) (2025 `size(PIE.Tu,2)`, `size(PIE.Tw,2)`),
  retagged 18 – equal to the prep verdicts.
- internal notes removed: none. Kept: `%    x_order_PDE = x_order(psize.n0+1:end);` (2025 commented code,
  `no` renamed `n0` inside it after 2025); three empty `% % ` lines of YP (no rule applies).
- open: 2025 entry continuation lines were changed by a global rename no -> n0 after 2025
  ("Terms format has n0w", "no0, nr0" for "noo, nro"); left as they are.

## 23 PIESIM/PIESIM_routines/PIESIM_Setup/PIESIM_options_check.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: ifexact is an option
  %   (opts.ifexact, default false), not a field of uinput; a PIE no longer
  %   needs the argument ndiff (opts.piesize removed; the piesize clause of
  %   6/16/2022 no longer applies); tf < dt with tf ~= 0 resets tf to 0 with
  %   a warning.
  ```
- replaced entries: MMP 11/19. The note marks the 2025 entry "piesize clause ... 6/16/2022" as no
  longer in the code (protocol: invalidated earlier entry noted at the end of the log).
- apply: 0/0/0. internal notes removed: none.

## 24 PIESIM/PIESIM_routines/PIESIM_Setup/rescalePIE.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: the maps PIE.misc.Tmap,
  %   Tumap, Twmap are translated (were PIE.misc.T0, Tu, Tw); help and
  %   comments revised.
  ```
- replaced entries: YP 1/6/2026, YP 4/7/2026, MMP 11/19 (2025 entry YP 6/26/2025, T0 -> Tmap
  reworded, stays). apply: 0/0/0. internal notes removed: none.

## 25 PIESIM/PIESIM_routines/PIESIM_Setup/rescalePIE_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: a given domain I is
  %   used (was replaced by the default whenever nargin < 3); more than 2
  %   inputs is an error; the maps PIE.misc.Tmap, Tumap, Twmap are translated
  %   (was PIE.misc.T0); help and comments revised.
  ```
- replaced entries: YP 1/6/2026, YP 4/7/2026, MMP 11/19 (2025 entry YP 6/26/2025 stays).
  apply: 0/0/0. internal notes removed: none. EOL i/crlf w/crlf unchanged.

## 26 PIESIM/PIESIM_routines/PIESIM_Setup/transl_1D.m
- class condense; TAG `MMP, 11/21/2026`
- note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 27 PIESIM/PIESIM_routines/PIESIM_Setup/transl_2D.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments
  %   revised.
  ```
- replaced entries: YP 4/7/2026 ("redefined rescaled operator domain": `T.I=I;` added 04/07/2026
  and removed by YP 05/02/2026, so no code change remains), MMP 11/19. apply: 0/0/0.
- open: YP kept in the TAG as a post-2025 committer although the net change of YP is nil; drop YP
  if the TAG should list only authors of surviving changes.

## 28 PIESIM/PIESIM_routines/PIESIM_Setup/transl_3PI.m
- class condense; TAG `MMP, 11/21/2026`
- note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 29 PIESIM/PIESIM_routines/PIESIM_Setup/transl_mod.m
- class condense; TAG `YP, MMP, 11/21/2026` (YP: removal of `T.I = I;`, 05/02/2026, no note)
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: T may be a 'struct' of
  %   its kernels as well as an 'opvar' (error message revised); the field
  %   T.I is no longer set to I; help and comments revised.
  ```
- replaced entries: MMP 11/04/2026 (2 lines), MMP 11/19.
- apply: dropped 0, kept 2 (was) (2025 opvar check), retagged 3. internal notes removed: none.

## 30 PIESIM/PIESIM_routines/PIESIM_Setup/vec_order2elem_order.m
- class condense; TAG `MMP, 11/21/2026`
- note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- replaced entries: MMP 11/19. apply: 0/0/0. internal notes removed: none.

## 31 PIESIM/PIESIM_routines/Postprocessing/PIESIM_plot_solution.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: solution fields ode, pde
  %   replaced by primary{1,2}, and observed, regulated are cells with the
  %   finite-dimensional outputs in {1}; psize fields n0, nr0, no0 (were no,
  %   nro, noo); grid is the array of grid points (was grid.phys);
  %   opts.ifexact replaces uinput.ifexact; printed values name the
  %   finite-dimensional outputs.
  ```
- replaced entries: YP 01/02/2026 (2 lines), YP 02/17/2026 (3), MMP 11/19 (2025 entries YP 6/16/2022
  and YP 06/02/2025, punctuation changed, stay). apply: 0/0/0. internal notes removed: none.

## 32 PIESIM/PIESIM_routines/Postprocessing/PIESIM_transform_to_solution.m
- class condense; TAG `VJ, YP, MMP, 11/21/2026` (VJ 12/29/2025: PIE.misc.Tu, Tw used when present;
  since rewritten by YP, but VJ is a post-2025 author)
- note:
  ```
  % VJ, YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: input grid removed;
  %   solution fields ode, pde replaced by primary{1,2}, and observed,
  %   regulated are cells {1,2} (finite-dimensional outputs; outputs of s in
  %   physical space); the input terms enter through Dop.Tucheb_2PDEstate,
  %   Twcheb_2PDEstate (or their map versions) with uinput.w, uinput.u cells
  %   of function handles of t; Dop.Tcheb_2PDEstate (or Tchebmap_2PDEstate)
  %   replaces Mcheb_nonsquare; psize fields n0, nr0, no0, nrx, nox; for a
  %   diverged (NaN) solution timedep.primary{2} is empty (was timepde.pde).
  ```
- replaced entries: YP 12/31/2025 (3 lines), YP 02/17/2026 (3), MMP 11/19, MMP 11/20/2026 (2)
  (2025 entry YP 6/16/2022, "not" -> "n0t" by the rename, stays).
- apply: dropped 0, kept 2 (was) (2025 code `if (ns>0)`, `solution.timepde.pde=[];`), retagged 2.
- internal notes removed: none.
- open (code, not fixed): the user message 'Warning: solution final time does n0t match user input
  final time' carries the rename typo "n0t" in a string literal (only this file in PIESIM/).

## 33 PIESIM/PIESIM_routines/Third_Party_Supporting_Files/fcgltran2d.m
- class condense (third-party, no license block added); TAG `YP, 11/21/2026`
- note (no note before; placed after `% Initial coding unknown`, the end of the log block, so that
  `help` does not print it):
  ```
  % YP, 11/21/2026: Changes since PIETOOLS 2025: spectral-to-nodal transform
  %   of N x M data with N ~= M corrected (the second transform keeps M
  %   values, was N).
  ```
- apply: 0/0/0. internal notes removed: none.

## 34 PIESIM/PIESIM_routines/Third_Party_Supporting_Files/fcgltran2dperm.m
- class condense (third-party, no license block added); TAG `YP, 11/21/2026`
- note: same text and place as 33. apply: 0/0/0.

## 35 PIESIM/PIESIM_routines/Third_Party_Supporting_Files/fcgltran4d.m
- class condense; TAG `MMP, 11/21/2026` (written by YP, not third-party; license kept)
- note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- replaced entries: MMP 11/19. apply: 0/0/0.

## 36 PIESIM/PIESIM_routines/Time_Integrators/PIESIM_int_bdf.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: new input gridall;
  %   uinput.w, wdot, u, udot are cells of function handles of t for separable
  %   inputs; a non-separable input (2D) enters through the coefficients of
  %   its spatial profile, recomputed at each step; Dop.Tcheb_inv replaces
  %   Dop.Mcheb_inv; psize.n0 replaces psize.no.
  ```
- replaced entries: YP 12/31/2025 (2 lines), YP 2/12/2026, MMP 11/19 (2025 entry YP 8_18_2024
  stays). apply: 0/0/0. internal notes removed: none.

## 37 PIESIM/PIESIM_routines/Time_Integrators/PIESIM_int_symbolic.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: the inputs w, wdot, u,
  %   udot are cells of function handles of t, evaluated at st in the
  %   convolution integrals; Dop.Tcheb_inv replaces Dop.Mcheb_inv; psize.no
  %   is not read; solcoeff.w, solcoeff.u are returned.
  ```
- replaced entries: MMP 11/19, MMP 11/20/2026 (2 lines, cited "YP 12/31/2025").
- apply: dropped 0, kept 5 (was) (2025 code), retagged 8 – equal to prep.
- comments rewritten: inline comment of `solcoeff.w/u=coeff.w/u` named the reading function
  ("read by PIESIM_transform_to_solution", release rule 6) -> "spatial coefficients of the inputs,
  for the reconstruction of the solution"; help line of uinput said "indexed as arrays and symbolic
  in st" (stale since the 11/20 fix) -> "cells of function handles of t".

## 38 PIESIM/PIESIM_routines/Time_Integrators/PIESIM_stability_check.m
- class condense; TAG `MMP, 11/21/2026`
- note: `% MMP, 11/21/2026: Changes since PIETOOLS 2025: help and comments revised.`
- replaced entries: MMP 11/19. apply: 0/0/0.

## 39 PIESIM/PIESIM_routines/Time_Integrators/PIESIM_time_integrate.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: new input gridall,
  %   passed to the BDF integrator, also when symbolic integration falls back
  %   to BDF2.
  ```
- replaced entries: YP 2/12/2026, MMP 11/19, MMP 11/20/2026 (2 lines).
- apply: dropped 0, kept 1 (was) (2025 call without gridall), retagged 1 (line longer than 76: marker
  one blank after the code, as cpass allows).

## 40 PIESIM/examples_pde_library_PIESIM_1D.m
- class condense; TAG `YP, MMP, 11/21/2026` (prep proposed DJ, VJ too: their lines `DJ, 12/19/2024:`
  and `VJ, 11/22/2025:` are 2025 entries reworded after 2025, so they stay; DJ, VJ have no commit
  in this file since 2025)
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: initial conditions in
  %   the single list uinput.ic (was uinput.ic.PDE, uinput.ic.ODE); a
  %   disturbance or control input may be a sym in st or a 2 x m data array
  %   [t; value] in a cell (example 2 uses data arrays; the alternative form
  %   is commented in the examples); the cylindrical examples 39-40
  %   (variable substitution) removed, 41-43 (weighted) renumbered 39-41,
  %   with uinput.weight only (uinput.weighted removed); example 1 on the
  %   domain [0.2,0.9] (was [-1,1]).
  ```
- replaced entries: YP 12/29/2025 (2 lines), YP 2/17/2026, MMP 11/19. apply: 0/0/0.
- comments rewritten: "accomplished intrinsically in PIESIM via weigthPIE routine" -> "... in PIESIM"
  (3 identical lines; release rule 6, no function names in comments).
- kept: the commented "Symbolic disturbance" / "Disturbance as a data array (double)" alternatives
  (new since 2025, unmarked) – deliberate user-facing alternatives of the input format (example 2
  runs the data-array form with the symbolic form commented), not development steps; and the
  commented terms-format, ic and batch blocks, which are 2025 commented code (fields renamed inside).
  Maintainer may prefer to drop the alternatives.

## Rerun with the fixed cpass.pl (coordinator, 11/21/2026)
- All prep runs above used a private copy of cpass.pl writing its diff files to the G1 scratch folder,
  so no prep diff was shared with another agent. After the fix, prep, apply and check were rerun with
  the shared tool on files 01-40: code diffs identical to the earlier ones (md5 of the CODE DIFF
  section), no marker other than the TAG found ("DJ (CC)" or no-blank markers: none in these files),
  apply idempotent (dropped 0; kept/retagged equal to the counts above), check PASS on 39 files and the
  2 known false positives on 18, codeident OK, EOL and numstat unchanged on all 40. No note changed.

## 41 PIESIM/examples_pde_library_PIESIM_2D.m
- class condense; TAG `YP, MMP, 11/21/2026` (prep proposed DJ too: `DJ, 01/04/2024:` is the 2025
  entry with ", update Example 19;" dropped after 2025; it stays; no DJ commit since 2025)
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: initial conditions in
  %   the single list uinput.ic (was uinput.ic.PDE, uinput.ic.x); 32 examples
  %   (was 19): new example 2 (inhomogeneous Dirichlet-Dirichlet), example 5
  %   (2025 example 2 with inhomogeneous Dirichlet-Neumann conditions), and
  %   examples 21-32 (boundary forcing, transport, vortex advection); the
  %   2025 examples 5-19 are 6-20.
  ```
- replaced entries: YP 2/15/2026, YP 2/17/2026, MMP 11/19. apply: 0/0/0. (prep, apply, check with the
  fixed tool; code diff identical to the earlier prep.)
- kept: commented example 9.5 (2025 commented code, was "case 8.5"), the commented PDE descriptions
  and the commented "Separable disturbances" alternative of example 26 (user-facing).

## 42 PIESIM/solver_PIESIM.m
- class condense; TAG `YP, MMP, 11/21/2026`
- note:
  ```
  % YP, MMP, 11/21/2026: Changes since PIETOOLS 2025: options are set before
  %   the example (opts.ifexact replaces uinput.ifexact); dim other than 1 or
  %   2 is an error; examples 1-41 in 1D (was 1-38) and 1-32 in 2D (was
  %   1-19); defaults example 41, tf = 0.1, dt = 1e-3 (were 1, 1.6, 0.02); a
  %   PIE is simulated without the argument ndiff; the PDE case also stores
  %   grid.
  ```
- replaced entries: YP 12/21/2025 (3 lines; its item "(2) pass options to examples" is not in the
  code – the example calls are unchanged since 2025), MMP 11/19. apply: 0/0/0 (fixed tool).
