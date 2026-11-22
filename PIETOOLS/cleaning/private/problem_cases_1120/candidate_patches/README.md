# Candidate example changes (problem-case run of 11/20/2026)

Diffs of the worktree copies (C:/wt1120b at 41c9210c) against 41c9210c, as run through the library path (libcheck):

- PDE_39.patch – PIETOOLS_PDE_Ex_2D_NS_Antonelli.m: form 'P', st.settings_2d = settings_2d_stock('light'),
  st.lpi.ineq.dD = 0 (the 'heavy' variant passed the library path as well).
- PDE_70.patch – PIETOOLS_PDE_Ex_2D_Wave_Acoustic_Boundary_Beale.m: the mixed position–strain form, st.lpi.ineq.dD = 0.
- PDE_69.patch – PIETOOLS_PDE_Ex_Wave_Heat_Transmission_2D_Batty.m: the mixed form ('69mix2'), st.lpi.ineq.dD = 0.
- PDE_44.patch, PDE_X3.patch – PIETOOLS_PDE_Ex_Wave_Eq_w_Boundary_Delay(_2D).m: the minimal form, the dual test
  (stability_dual with form 'P'), st.epneg = 1e-1 (44) or 1e-2 (X3).
- all_examples_worktree_41c9210c.patch – every example edit of that worktree, including the TRADE3 statements and
  the checks of 36 and 41.

The marker "(worktree test)" on each line is the worktree's; an edit of the release files carries its own header
entry and markers. Not applied: the release comment pass edits these files, and the formulation changes (44, X3, 69,
70) are the maintainer's call.
