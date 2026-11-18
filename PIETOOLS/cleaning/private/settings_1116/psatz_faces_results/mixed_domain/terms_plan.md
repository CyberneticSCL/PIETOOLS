# Term experiments on the mixed suite and the row tiers: plan (2026-11-16, 05:10)

Maintainer, 11/16 05:05: "Start the term experiments on the mixed suite. Run autonomously as
yesterday for 12 hours. At the end, create a pdf report detailing the construction, the
proposed knobs as related to the construction, and the proposed settings for each knob in the
new proposed tiers. The new tiers should approximately replicate the ability of existing tiers
to handle single domains without big significant increases in complexity, but should also
perform comparably well on mixed domains."

Window: 05:10 to 17:10. Rules as yesterday: scratch shadows only, nothing in the repository
changed, every batch logged, decisions in `../DECISIONS.md` (D37 on), results in
`terms_results.md` here, the report `sopvar/settings_report_2026-11-16.tex` (+ pdf) at the end.
MOSEK at 8 threads, one MATLAB process at a time, batches of about 20-30 minutes with their
own logs under the session scratchpad `psatz_faces/logs/t*.log`.

## Instruments (Phase 0)

- `pf_row_storage.m`: the proposal's reader. A row per entry (number of variables of the
  space): len_inc d1, len_com d2, len_sos d3 (rows a, b, c), patterns (uniform, separable,
  all, or explicit), terms and offsets (storage: none). Basis operators: patterns over
  {a, b, c} = codes {1, 4, 2} of poscopvar_direct; per operator mult_k = length(alpha_k) - 1
  (0, d2_k - 1, d1_k - 1), int_k = d3_{alpha_k,k} - 1; no total cap; R^q: identity at its
  multiplier length. Built by POSCOPVAR with explicit DEG and include.
- `route_row/exec_tools_sop.m`: the executives' storage in 1-D, 2-D and N-D from the row
  (appdata PF_ROW); slack, caps and margins unchanged.
- `pf_qdirect.m` with term options on both inequalities; `pf_poinc.m` with term options.
- Term vocabulary (psatz_rows over the registry): none; faces (scalar offset, 'auto');
  faces at offset 0; faces per direction; faces of one direction only; box (Putinar, per
  direction); box of one direction; corners.

## Phase 1: the terms (06:00-08:30)

T1 the lift of the PDE norm, gam I - T'T, reach by term set (S3, S8, S17, S18, S26).
T2 the direct test at 0.5 and 0.9 of the margin by term set, Q at length 2 (S3, S8, S17, S18,
   S26, S4, S6); the executive 'P' form with the row storage where K is present (S9, S24).
T3 offsets scalar against per direction, at w = 1 and w = 2, on T1 and T2.
T4 the 1-D space's own box factor as a row, alone and with the faces.

## Phase 2: the row tiers against the baselines (08:30-11:30)

Candidate rows R-light, R-heavy (and R-veryheavy where the stock tier is run): single-domain
baselines (1-D: program-level comparison with the stock storage on library rows, Prop. span;
2-D: rows 33-35, 37, 39, 40, 45-48 and the fin2 ladder; 3-D: heatNd light stability and gain),
then the mixed suite (exec layer, 41 checks) and S3/S8/S18 at 0.9, S17 by the direct test.

## Phase 3: adjustment (11:30-14:00), Phase 4: the report (14:00-17:00)

## Additions (maintainer, 11/16 06:35)

- **Per-direction caps as the knob.** In place of any total-degree or total-length cap (the
  N-D `joint` cap, the 2-D cap tables), a cap per direction on the degree of the kernel
  monomials in that direction, for every block and every space (the triangular grid
  theta_k + s_k <= c_k that `opts.cap` of lpi_ineq_sop gives the slack on a full-registry
  space, extended to the storage and to spaces of fewer variables). Measurement: no cap
  against c_k against the 2-D tables at equal verdict on the sensitive baselines (fin2, PDE
  48, the library rows), in size and time.
- **Revision of poscopvar in scope.** Where the row needs it (per-direction caps on a
  lower-dimensional space; the pattern sets with incomplete moments; the non-owned
  direction rule), the change to poscopvar / poscopvar_direct is part of the agenda, in a
  scratch shadow first and proposed as a diff.
- **poscopvar_lift back on the path.** poscopvar is the efficient implementation of
  poscopvar_lift (`orphaned_file/private/poscopvar_lift.m`, the explicit A^* M A of the
  construction with D, w per direction); it is to be moved back into the path as an
  analysis tool (working tree move, no commit).
