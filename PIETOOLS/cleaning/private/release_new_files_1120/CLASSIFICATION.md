# New files of the 2026 release: renames, replacements, new functionality

20 November 2026, HEAD 6c2cd695 (ndopvar). Full table: `new_files_classified.tsv` (this folder).

## Scope

The list is the census of the release comment pass of the Repo cleanup session
(`cleaning/private/release_comment_pass/lists/`, taken 19 November 2026, 22:46): the files marked
new in the six group files, and the three `*_newfix.txt` lists of class methods the first census
passed as existing. That gives 567 files. The same rule applied at HEAD adds 26 files the census does
not hold: 24 it deferred because other sessions had uncommitted edits in them (21 library examples,
`DEMO15_trace_inequality_2D`, `PIETOOLS_auto_execute_sop`), and 2 PIESOS files added after it
(`get_tensorvars`, `@quadPoly/times`, Declan, 19 November). Total: 593 files.

A file is new when its key is absent from the 2025 public release `origin/PIETOOLS-2025`
(fb7f1209, 29 December 2025): the key is `@class/method` inside a class folder, the basename
elsewhere. This is the rule of the census after its newfix correction. The census covers the release
folders only (no SOSTOOLS400, testing_codes, cleaning, PIETOOLS_demos/cuadmm, lowrank_2d_stability).

## Classes

- **Rename**: the code descends from a 2025 file that did the same thing. Evidence: git rename
  detection between the 2025 tree and HEAD (similarity in the table), or a block of a 2025 file moved
  into its own file.
- **Replacement**: no code descends from 2025, but the file does what a named 2025 file did. Most
  replace a 2025 file for the operator classes that took over from `opvar`/`dopvar`/`opvar2d`/`dopvar2d`
  (`copvar`, `cdopvar`, `sopvar`, `sdopvar`), or for container PIEs. Several do more than their 2025
  counterpart: any number of spatial variables, spaces of mixed dimension, or one file in place of
  several (the executives and the settings tiers).
- **New**: no 2025 file does the same thing.

Level: *internal* is a file under `private/`, `misc/`, `internal/`, `utility_functions/`,
`PIESIM_routines/`, `converters/legacy/`, or third-party; *user* is every other file.

## Counts

| folder | Rename | Replacement | New, user | New, internal | total |
|---|---|---|---|---|---|
| PIESIM | 14 | 0 | 0 | 3 | 17 |
| PIESOS | 0 | 0 | 206 | 31 | 237 |
| PIETOOLS_demos | 0 | 0 | 7 | 0 | 7 |
| PIETOOLS_examples | 0 | 0 | 36 | 0 | 36 |
| converters | 1 | 1 | 14 | 13 | 29 |
| executives_sopvar | 0 | 13 | 1 | 14 | 28 |
| lpi_programming_sopvar | 0 | 9 | 1 | 19 | 29 |
| parser | 0 | 0 | 13 | 10 | 23 |
| settings_sopvar | 0 | 8 | 1 | 0 | 9 |
| sopvar | 0 | 82 | 48 | 48 | 178 |
| **total** | **15** | **113** | **327** | **138** | **593** |

So 128 files (15 renames, 113 replacements) do what a 2025 file did, and 465 are new: 327 at the
user level, 138 internal. PIESOS accounts for 237 of the new files; it is not in the 2025 release.

## Renames (15)

| new file | 2025 file | evidence |
|---|---|---|
| `PIESIM_1Dto1D2Mat_cheb_2D.m` | `PIESIM_PI2Mat_cheb_opint_discretize_1to1_line.m` | git R061 |
| `PIESIM_1Dto2D2Mat_cheb_2D.m` | `PIESIM_1Dto2D2Mat_cheb.m` | git R047 |
| `PIESIM_2Dto0D2Mat_cheb_2D.m` | `PIESIM_PI2Mat_cheb_opint_discretize_2to0.m` | git R069 |
| `PIESIM_2Dto1D2Mat_cheb_2D.m` | `PIESIM_2Dto1D2Mat_cheb.m` | git R049 |
| `PIESIM_3PI_op_discretize_2D.m` | `PIESIM_PI2Mat_cheb_opmixed_discretize_1to2.m` | git R043 |
| `PIESIM_3PI_op_line_int_discretize_2D.m` | `PIESIM_PI2Mat_cheb_opmixed_discretize_2to1.m` | git R041 |
| `PIESIM_3PI_opint_discretize.m` | `PIESIM_3PI2Mat_cheb_opint_discretize.m` | git R088 |
| `PIESIM_3PI_opmult_discretize.m` | `PIESIM_PI2Mat_cheb_opmult_discretize.m` | git R057 |
| `PIESIM_9PI_opint_area_discretize_2D.m` | `PIESIM_9PI2Mat_cheb_opint_discretize_area.m` | git R064 |
| `PIESIM_9PI_opint_line_discretize_2D.m` | `PIESIM_9PI2Mat_cheb_opint_discretize_line.m` | git R058 |
| `PIESIM_9PI_opmult_discretize_2D.m` | `PIESIM_PI2Mat_cheb_opmult_discretize_2D.m` | git R057 |
| `PIESIM_PI2Mat_opint_cheb.m` | `PIESIM_PI2Mat_cheb_opint_discretize.m` | git R091 |
| `PIESIM_PI_discretize.m` | `PIESIM_PIop_discretize.m` | git R083 |
| `PIESIM_disturbance_check.m` | `PIESIM_input_check.m` | split out: the disturbance check and split of the 2025 input check, which now calls it |
| `convert_PIETOOLS_PDE_2D.m` | `convert_PIETOOLS_PDE.m` | the 2025 converter (git C074), kept as the 2-D path; convert_PIETOOLS_PDE.m is now the N-D converter |

## Replacements (113)

### Functions (35): 2025 name and new name

| new file | 2025 file(s) | what it does |
|---|---|---|
| `PIETOOLS_H2_control_sop.m` | PIETOOLS_H2_control.m | H2-optimal state feedback; same executive for container PIEs |
| `PIETOOLS_H2_estimator_sop.m` | PIETOOLS_H2_estimator.m | H2-optimal estimator; same executive for container PIEs |
| `PIETOOLS_H2_norm_sop.m` | PIETOOLS_H2_norm_c.m, PIETOOLS_H2_norm_c_coercive.m, PIETOOLS_H2_norm_o.m, PIETOOLS_H2_norm_o_coercive.m, PIETOOLS_H2_norm_2D_c.m, PIETOOLS_H2_norm_2D_c_non_coercive.m, PIETOOLS_H2_norm_2D_o.m, PIETOOLS_H2_norm_2D_o_non_coercive.m | one executive in place of eight (form argument; any dimension) |
| `PIETOOLS_Hinf_control_sop.m` | PIETOOLS_Hinf_control.m | H-infinity state feedback; same executive for container PIEs |
| `PIETOOLS_Hinf_estimator_sop.m` | PIETOOLS_Hinf_estimator.m, PIETOOLS_Hinf_estimator_2D.m | one executive in place of the 1-D and 2-D estimators |
| `PIETOOLS_Hinf_gain_sop.m` | PIETOOLS_Hinf_gain.m, PIETOOLS_Hinf_gain_coercive.m, PIETOOLS_Hinf_gain_dual.m, PIETOOLS_Hinf_gain_dual_coercive.m, PIETOOLS_Hinf_gain_2D.m, PIETOOLS_Hinf_gain_2D_non_coercive.m, PIETOOLS_Hinf_gain_dual_2D.m | one executive in place of seven (form argument; any dimension) |
| `PIETOOLS_auto_execute_sop.m` | PIETOOLS_auto_execute.m | runs the executives that workspace flags select (GUI/script path) |
| `PIETOOLS_stability_sop.m` | PIETOOLS_PDEstability.m, PIETOOLS_PDEstability_dual.m, PIETOOLS_PIE2PDEstability.m, PIETOOLS_PIE2PDEstability_dual.m, PIETOOLS_stability_2D.m, PIETOOLS_stability_dual_2D.m | one executive in place of six (form argument; any dimension) |
| `PIETOOLS_well_posedness_sop.m` | PIETOOLS_well_posedness.m | well-posedness test; same executive for container PIEs |
| `apply_sopvar.m` | apply_opvar.m | applies an operator to a function |
| `degbalance_core.m` | @opvar/degbalance.m, @opvar2d/degbalance.m | shared core of the degbalance methods |
| `eq_copvar.m` | @opvar/eq.m, @opvar2d/eq.m | shared core of the eq methods of the containers |
| `getController_sop.m` | getController.m | controller K from P, Z |
| `getObserver_sop.m` | getObserver.m, getObserver_2D.m | observer gain L from P, Z |
| `get_lpivar_degs_sop.m` | get_lpivar_degs.m | degrees of the free operators |
| `getsol_lpivar_sop.m` | getsol_lpivar.m, getsol_lpivar_2d.m | solved value of a decision operator |
| `lpi_eq_cdopvar.m` | lpi_eq.m | container part of lpi_eq_sop |
| `lpi_eq_sdopvar.m` | lpi_eq.m | block part of lpi_eq_sop |
| `lpi_eq_sop.m` | lpi_eq.m, lpi_eq_2d.m | operator equality constraint |
| `lpi_ineq_sop.m` | lpi_ineq.m, lpi_ineq_2d.m | operator inequality constraint |
| `lpigetsol_sop.m` | lpigetsol.m | solved value of a decision variable |
| `lpiscript_sop.m` | lpiscript.m | predefined LPIs by name |
| `lpisettings_sop.m` | lpisettings.m | settings by tier name |
| `lpivar_cdopvar.m` | lpivar.m, lpivar_2d.m | free decision operator |
| `mat2copvar_sop.m` | mat2opvar.m | matrix to (container) multiplier operator |
| `mat2sopvar.m` | mat2opvar.m | matrix to (block) multiplier operator |
| `poscopvar.m` | poslpivar.m, poslpivar_2d.m | positive decision operator |
| `poslpivar_sop.m` | poslpivar.m | positive operator Z'(gQ)Z as the executives declare it |
| `settings_2d_stock.m` | settings_PIETOOLS_*_2D.m | the 2-D settings of the stock tiers, by level |
| `settings_sopvar_custom.m` | settings_PIETOOLS_custom.m, settings_PIETOOLS_custom_2D.m | tier file (1-D and 2-D files merged) |
| `settings_sopvar_extreme.m` | settings_PIETOOLS_extreme.m, settings_PIETOOLS_extreme_2D.m | tier file (1-D and 2-D files merged) |
| `settings_sopvar_heavy.m` | settings_PIETOOLS_heavy.m, settings_PIETOOLS_heavy_2D.m | tier file (1-D and 2-D files merged) |
| `settings_sopvar_light.m` | settings_PIETOOLS_light.m, settings_PIETOOLS_light_2D.m | tier file (1-D and 2-D files merged) |
| `settings_sopvar_stripped.m` | settings_PIETOOLS_stripped.m, settings_PIETOOLS_stripped_2D.m | tier file (1-D and 2-D files merged) |
| `settings_sopvar_veryheavy.m` | settings_PIETOOLS_veryheavy.m, settings_PIETOOLS_veryheavy_2D.m | tier file (1-D and 2-D files merged) |

### Class files and methods (78)

A method counts as a replacement when the 2025 class of the same role has a method of the same name:
`@opvar`/`@opvar2d` for `copvar` and `sopvar`, `@dopvar`/`@dopvar2d` for `cdopvar` and `sdopvar`.
The class files `copvar.m` and `cdopvar.m` replace `opvar`/`opvar2d` and `dopvar`/`dopvar2d`.

| class | files | methods |
|---|---|---|
| `@cdopvar` | 17 | blkdiag, cdopvar, ctranspose, display, end, eq, horzcat, minus, mtimes, numArgumentsFromSubscript, op_slice, plus, size, subsasgn, subsref, uminus, vertcat |
| `@copvar` | 20 | blkdiag, copvar, ctranspose, diff, display, end, eq, horzcat, inv, minus, mtimes, numArgumentsFromSubscript, op_slice, plus, size, subs, subsasgn, subsref, uminus, vertcat |
| `@sdopvar` | 17 | blkdiag, ctranspose, degbalance, display, end, eq, horzcat, minus, mtimes, numArgumentsFromSubscript, op_slice, plus, size, subsasgn, subsref, uminus, vertcat |
| `@sopvar` | 24 | blkdiag, clean_op, ctranspose, degbalance, diff, display, end, eq, horzcat, inv, inv_1D, minus, mtimes, numArgumentsFromSubscript, op_slice, plus, size, sopvar2sdopvar, subs, subsasgn, subsref, uminus, uplus, vertcat |

## New functionality, user level (327)

| folder | files | files (names) |
|---|---|---|
| PIESOS | 40 | Liediff, DMB, SOS_DP, filter_Pcell, kernel_degree, polyopvar_equal, weighted_sobolev_norm, PIESOS_DEMO1_global_stability, PIESOS_DEMO2_local_stability, PIESOS_DEMO3_energy_bound, examples_PDE_library_PIESOS, dmonomials, PIESOS_settings, PIESOS_stability_global, PIESOS_stability_local, PIESOS_sublevel, innerprod, monopmat, piesos_LF, piesos_eq, piesos_getsol, piesos_ineq, piesos_intopvar, piesos_lpivar, piesos_polyvar, piesos_poslpivar, piesos_program, piesos_quadopvar, piesos_quadvar, piesos_setobj, piesos_solve, piesos_sosvar, constant2polyopvar, copvar2tensopmat, mat2tensopmat, mat2tensopvar, pdevar2polyopvar, sopvar2tensopvar, tensopmat2copvar, tensopvar2polynomial |
| PIESOS/polyopvar/@intop | 31 | apply_op, blkdiag, canonicalize_intmat, clean_op, combine_terms, common_intmat, ctranspose, diff, display, get_display, horzcat, intop, is_id_op, isvalid, kernels, minus, mtimes, mtimes_tensopvar, numArgumentsFromSubscript, numel, otimes, plus, reorder_states, reorder_vars, size, subsref, times, transpose, uminus, unique_intmat, vertcat |
| PIESOS/polyopvar/@polyopvar | 35 | apply_poly, clean_poly, combine_factors, combine_terms, common_basis, common_outputs, common_states, common_vars, ctranspose, diff, display, get_display, horzcat, initialize, ispvar, lift_output, minus, mpower, mtimes, num_factors, opower, otimes, plus, polyopvar, power, remove_unused_states, repelem, repmat, size, subs, subsasgn, subsref, times, uminus, vertcat |
| PIESOS/polyopvar/@quadPoly | 38 | Sym2quadPoly, blkdiag, combine, common_vars, ctranspose, decomposeLR, diff, disp, end, eq, hessian, horzcat, int, isempty, jacobian, minus, mtimes, numArgumentsFromSubscript, numel, plus, polynomial2quadPoly, quadPoly, quadPoly2polynomial, randquadPoly, repmat, reshape, set_vars, size, subs, subsasgn_1, subsref, trace, transpose, uminus, uplus, var_swap, vertcat, times |
| PIESOS/polyopvar/@tensopmat | 28 | apply_op, blkdiag, clean_op, common_vars, ctranspose, diff, display, get_display, horzcat, initialize, is_id_op, isvalid, lift_output, minus, mtimes, opower, otimes, plus, repelem, repmat, size, subsasgn, subsref, tensopmat, times, transpose, uminus, vertcat |
| PIESOS/polyopvar/@tensopvar | 34 | apply_op, blkdiag, clean_op, combine_factors, common_vars, ctranspose, diff, display, eq, get_display, horzcat, initialize, int, is_id_op, isvalid, kron, lift_output, minus, mtimes, num_factors, opower, otimes, plus, reorder_vars, repelem, repmat, size, subsasgn, subsref, tensopvar, times, transpose, uminus, vertcat |
| PIETOOLS_demos | 7 | DEMO10_PBC, DEMO11_poincare_inequality_2D, DEMO12_2DControl, DEMO13_Hinf_OFcontrol, DEMO14_trace_inequality, heatNd_apply, DEMO15_trace_inequality_2D |
| PIETOOLS_examples | 36 | l2gain_23, l2gain_26, stability_13_2, stability_13_3, stability_16_1, stability_16_3, stability_16_4, stability_17, stability_18, stability_19, stability_36, stability_39, stability_42, examples_FI_library_PIETOOLS, PIETOOLS_PDE_Ex_2D_Plate_Eq, PIETOOLS_PDE_Ex_2D_Wave_Eq_Damped, PIETOOLS_PDE_Ex_3D_Advection_Diffusion_Eq, PIETOOLS_PDE_Ex_3D_Coupled_Heat_Eq, PIETOOLS_PDE_Ex_3D_Reaction_Diffusion_Eq, PIETOOLS_PDE_Ex_2D_Delayed_Heat_Eq_Fridman, PIETOOLS_PDE_Ex_2D_Heat_Eq_Exact_Gain_Lidstrom, PIETOOLS_PDE_Ex_2D_Linear_Turing_Pair, PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale, PIETOOLS_PDE_Ex_Age_Structured_Diffusion_Kang, PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Fridman, PIETOOLS_PDE_Ex_Damped_Wave_Interior_Delay_Riemann_Fridman, PIETOOLS_PDE_Ex_Elastic_Thermoelastic_Transmission_Marzocchi, PIETOOLS_PDE_Ex_Heat_Delayed_Robin_Boundary_Nicaise, PIETOOLS_PDE_Ex_Heat_Edge_Cascade_Variable_Coeff, PIETOOLS_PDE_Ex_Heat_Laplace_Beltrami_Edge, PIETOOLS_PDE_Ex_Road_Field_Berestycki, PIETOOLS_PDE_Ex_Transport_Edge_Diffusion, PIETOOLS_PDE_Ex_Two_Layer_Heat_Rectangles, PIETOOLS_PDE_Ex_Wave_Delayed_Boundary_Damping_Xu, PIETOOLS_PDE_Ex_Wave_Heat_Transmission_1D_Batty, PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty |
| converters | 14 | convert_PIETOOLS_BVP, dopvar2sdopvar, dpvar2sdvar, dpvars2sdvars, npie2pie_struct, op2copvar_sop, opvar2copvar, opvar2d2sopvar, opvar2sopvar, pie2copvar, pie2npie_struct, pie2piesim, sdvar2dpvar, print_convert_summary |
| executives_sopvar | 1 | PIETOOLS_Hinf_OFcontrol_sop |
| lpi_programming_sopvar | 1 | lpisetsettings |
| parser/@npie_struct | 8 | display, eq, expand_tderivatives, horzcat, initialize, npie_struct, subsasgn, vertcat |
| parser/@pde_struct | 5 | mpower, power, static_constraints, times, times_pdes |
| settings_sopvar | 1 | lpisettings_stock_sop |
| sopvar | 5 | eye_copvar_sop, int2copvar, int2sopvar, pi4_copvar_sop, zeros_copvar_sop |
| sopvar/@cdopvar | 11 | cat, get_display, getfield_opvar1d, getfield_opvar2d, is_opvar1d, is_opvar2d, plus_batch, setfield_opvar1d, setfield_opvar2d, transpose, verify |
| sopvar/@copvar | 13 | cat, get_display, getfield_opvar1d, getfield_opvar2d, int, is_id_op, is_opvar1d, is_opvar2d, plus_batch, setfield_opvar1d, setfield_opvar2d, transpose, verify |
| sopvar/@sdopvar | 8 | cat, get_display, plus_batch, repelem, repmat, sdopvar, setdvars, transpose |
| sopvar/@sopvar | 11 | cat, get_display, int, is_id_op, kron, lift_output, repelem, repmat, sopvar, transpose, unkron |

By subject:

- **PIESOS (206 user, 31 internal)**: SOS programming for nonlinear PIEs: the classes `polyopvar`,
  `tensopvar`, `tensopmat`, `intop`, `quadPoly`; the program functions `piesos_*` (program, variables,
  constraints, objective, solve, solution); the executives `PIESOS_stability_global`,
  `PIESOS_stability_local`, `PIESOS_sublevel`; three demos and a PDE library.
- **Operator classes**: the block classes `sopvar`, `sdopvar` (one PI operator in any number of
  variables); container methods with no 2025 method: `cat`, `transpose`, `int`, `kron`, `repmat`,
  `repelem`, `plus_batch`, `verify`, `is_id_op`, `lift_output`, `unkron`, `setdvars`, and the
  opvar-style field access `getfield_opvar1d/2d`, `setfield_opvar1d/2d`, `is_opvar1d/2d`; the
  constructors `eye_copvar_sop`, `zeros_copvar_sop`, `int2copvar`, `int2sopvar`, `pi4_copvar_sop`
  (in 2025 an identity, a zero or a 4-PI operator was built by assigning the fields of an `opvar`).
- **Executives**: `PIETOOLS_Hinf_OFcontrol_sop` (H-infinity output-feedback control).
- **LPI programming and settings**: `lpisetsettings` (settings stored on the program),
  `lpisettings_stock_sop` (the stock slack for baselines).
- **Parser**: the class `npie_struct` (nonlinear PIEs); products and powers of PDE terms
  (`times`, `times_pdes`, `power`, `mpower`) and `static_constraints` for `pde_struct`.
- **Converters**: the N-D converter for boundary-value problems (`convert_PIETOOLS_BVP`); conversions
  between the 2025 classes and the new ones (`opvar2copvar`, `opvar2sopvar`, `opvar2d2sopvar`,
  `dopvar2sdopvar`, `dpvar2sdvar`, `dpvars2sdvars`, `sdvar2dpvar`, `op2copvar_sop`, `pie2copvar`);
  nonlinear PIEs (`pie2npie_struct`, `npie2pie_struct`); `pie2piesim`; `print_convert_summary`.
- **Demos (7)**: periodic boundary conditions (DEMO10), 2-D Poincare inequality (DEMO11), 2-D state
  feedback (DEMO12), H-infinity output feedback (DEMO13), trace inequalities in 1-D and 2-D (DEMO14,
  DEMO15), `sopvar_demos/heatNd_apply`.
- **Examples (36)**: the functional-inequality library `examples_FI_library_PIETOOLS`; 13 command-line
  examples; 2-D (plate, damped wave), 3-D (3) and mixed-domain (16) library examples.

## New, internal (138)

| folder | files | files (names) |
|---|---|---|
| PIESIM/PIESIM_routines | 3 | PIESIM_Mat_Truncate, weightPIE, clencurt |
| PIESOS/internal | 9 | clean_dpvar, coeff_error, common_vars_quadPolyCell, find_common_idcs, is_eye, otimes_params, quad_expand, set_vars_quadPolyCell, get_tensorvars |
| PIESOS/polyopvar | 22 | Zdims, Zmap, Zpermmap, display_term, enumerate_valid_orders, get_int_order, kernels2params, lmult_params, params2kernel, prune_basis, remap_cols, remap_dvars, sync_basis, factor_key, liftIndex, tensExps, unionBasis, is_native_op, merge_out_spaces, native2op, native_mask, native_vars |
| converters/internal | 10 | compose_Tmaps, compute_Tmap_1D, compute_canonical_BCs, decouple_BVPs, eliminate_BVP_sources, fund_inputs, lift_BVP_sources, split_BVPs, compute_PIETOOLS_PDE_periodic_map, convert_PDE_term |
| converters/legacy | 3 | compute_PIETOOLS_PDE_periodic_map_2D, convert_PIETOOLS_PDE_nonlinear_2D, pie_pde_maps_2D |
| executives_sopvar/utility_functions | 14 | exec_info_sop, exec_kmult_sop, exec_pie_sop, exec_sop_settings, exec_tools_sop, lpi_solve_loop_sop, pie_disc_sop, pie_witness_sop, poslpivar_settings_2d_sop, poslpivar_settings_sop, qp_copvar, settings2possopvar, settings_2d_sop, trace_rn_sop |
| lpi_programming_sopvar/misc | 19 | copquadvar, eq_opts_sopvar, get_lift_degs, lift_copvar, lpi_classify_sop, lpi_shape_sop, lpi_soseq, lpigetdual_sop, poscopvar_direct, poscopvar_lift, posmult_cdopvar, build_exponent_grid, collect_eq_rows, dvar_rows, expand_full, impose_eq_rows, merge_monomial_product, process_degrees_one, unpack_sheets |
| parser/@npie_struct | 3 | pie_compat_tab, pie_compat_vars, pie_space_canon |
| parser/@pie_struct | 4 | pie_compat_tab, pie_compat_vars, pie_rename_dummy, pie_space_canon |
| parser/internal | 3 | print_terms, use_unicode, varname2display |
| sopvar/@cdopvar | 2 | merge_dvar_lists, merge_dvar_pair |
| sopvar/@sdopvar | 9 | ChangeDecVar, apply_basis_map, decode_params, display_term, leftShiftMonomials_SD, leftShiftMonomials_SS, lrmultiply, params2kernel, sync_basis |
| sopvar/@sopvar | 5 | UnionBasisMonomials, decode_params, display_term, leftShiftMonomials_SS, params2kernel |
| sopvar/misc | 32 | canonical_adjoint_map, canonical_var_order, canonicalize_multiplier, is_canonical_multiplier, monomial_gather, prune_zero_monomials, cat_copvar_grid, chebfit_mono, compose_sopvar, cell_of_gamma, check_reserved_names, copvar_space_list, gamma_of_cell, kron_split, kron_strides, monomial_position, monomial_position_map, multiindex_grid, parse_copvar_spaces, derive_copvar_meta, dpvar_op_copvar, gk_grid, int_semisep, mat2copvar_grid, merge_copvar_registry, on_registry_sop, parse_diff_inputs, parse_int_inputs, pvar2varname, spIkron, split_copvar, verify_copvar_meta |

## Judgement calls (each can be overridden in the table)

1. `sopvar.m` and `sdopvar.m` are New (block classes; 2025 had no single-block N-D operator), while
   `copvar.m` and `cdopvar.m` are Replacements: they took over the role of `opvar`/`dopvar` and
   `opvar2d`/`dopvar2d`.
2. Internal cores that carry out a 2025 operation are Replacements: `degbalance_core`, `eq_copvar`,
   `@sopvar/inv_1D` (2025 `inv_opvar`), `@sopvar/clean_op` (2025 `@opvar/clean`), `lpi_eq_cdopvar`,
   `lpi_eq_sdopvar`, `poslpivar_sop` (the positive operator of the executives, 2025 `poslpivar`).
3. `PIETOOLS_stability_sop`, `PIETOOLS_Hinf_gain_sop`, `PIETOOLS_H2_norm_sop`,
   `PIETOOLS_Hinf_estimator_sop` and the six tier files each replace several 2025 files (listed).
4. `convert_PIETOOLS_PDE_nonlinear_2D` is New although git pairs it with the 2025
   `convert_PIETOOLS_PDE.m` (58 %): nonlinear conversion did not exist in 2025.
5. `stability_13_2` is New although it declares the beam of the 2025 `stability_13`, written as one
   equation second order in time.
6. The `npie_struct` methods are New (new class) although `pie_struct` has methods of the same names.
7. Git similarities of 40 – 65 % between short methods and 2025 files such as `@opvar/uminus.m` come
   from the shared license header and are not counted as evidence.
8. `PIESIM_disturbance_check` is a Rename: the 2025 `PIESIM_input_check` held this block and now calls
   the new file. `PIESIM_Mat_Truncate` is New: no 2025 file contains its operation.

## Table columns (`new_files_classified.tsv`)

`class` (Rename, Replacement, New), `level` (user, internal), `group` (census group),
`path` (under PIETOOLS/), `counterpart_2025`, `purpose` (help line or note), `list` (census, newfix,
deferred, after census).
