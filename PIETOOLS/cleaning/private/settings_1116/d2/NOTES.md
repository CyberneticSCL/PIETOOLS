# D2 at 0.5 tau* (library case 56), light 'Q', 3-D term entry 'faces1': sizes (MMP, 11/18/2026 05:53)

Size only (d2_size2.m with the LB_SIZE_ONLY gate; d2_size2.log). cap = sop.cap_delta(3) on every
slack basis (theta + s <= Jmax(T) + cap); cap_inc = lpi.ineq.cap_inc 'auto' = Dmax - 1 on the
lower and upper basis operators only (Shivakumar, Lemma 1.2: lossless for a purely integral target
with any Gram factorization; D2 has no ODE, its targets are purely integral).

| cap | cap_inc | variables | m | positivity slack | negativity slack |
|---|---|---|---|---|---|
| 1 | none (control, D81 rung) | 2,721,692 | 82,551 | 658,748 | 2,062,944 |
| 1 | auto | 1,394,588 | 63,066 | 658,748 | 735,840 |
| 0 | auto | 475,404 | 33,025 | 157,308 | 318,096 (cap 0 alone gives the same: cap_inc does not bind) |
| 2 | auto | 2,027,420 | 84,890 | 1,291,580 | 735,840 |
| none | auto | 2,134,428 | 90,777 | 1,398,588 | 735,840 (uncapped without cap_inc: 6,813,852 / 163,847) |

Reading: with cap_inc the negativity slack is the same 735,840 at every cap, so the rungs differ
only by the positivity slack; the uncapped program with the triangular cap (2,134,428 / 90,777) is
the principled rung: it keeps every certificate of the uncapped program (6.8M / 163,847, out of
memory on Sol at 200 GB) if the bound applies, at 3.2x fewer variables and 1.8x fewer rows. The
cap-1 control ran 4 MOSEK iterations on the workstation at about 190 s each with a sparse
factorization (d2_size2_cap1_mosek_partial.log), so the dense 23-25 bytes per m^2 model (D90) does
not govern these programs; a workstation solve of the 90,777-row rung is plausible in a few hours.
Scheduling is the Sister session's (study A3/A4).

## Memory of the cap-1 control (the accidental solve, 04:38-05:38; d2_size2_cap1_mosek_partial.log)
Not sampled. MOSEK's own statistics: 82,551 constraints, 74,378 after presolve, 14 matrix variables
(scalarized 1,362,168), 24 threads; factor nonzeros 1.40e9 before and 1.41e9 after factorization
(about 11.3 GB at 8 bytes; a dense factor at 74,378 rows would be 2.77e9), dense dim 3, 4.39e13
flops, setup 310 s, about 190 s per iteration. At launch 41 GiB were free with one other MATLAB
(6.6 GB) running, and the run never failed on memory, so its peak was below about 45 GB. The Sister's
uncapped + triangular-cap rung (m 90,777) was stopped by its 56 GB watchdog during the MOSEK setup
(17 -> 57 GB in two minutes), i.e. that factor is nearly dense (90,777^2/2 x 8 bytes = 33 GB for the
factor alone). The cap 1 + triangular cap rung (m 63,066; a dense factor would be 15.9 GB) should fit.

## The storage form 'P' at 0.5 tau* (d2_size3.m, size only, 07:12; the tags did not print, order as listed)
| form, tier | cap_inc | variables | m | note |
|---|---|---|---|---|
| 'P', light | none | 3,893,508 | 95,216 | slack like R: certificate 3,885,408, Gram blocks [1872 252 x6]; no cheaper than 'Q' uncapped + triangular cap (2,134,428 / 90,777) |
| 'P', light | auto | 9,108 | 12,826 | certificate 1,008 (Gram blocks [12 x7]): the cap Dmax - 1 is taken from the like-operator R, not from the target T'PA + A'PT, so Lemma 1.2 gives no guarantee here; a certificate would still count (one-minute solve), an infeasibility would not |
| 'P', heavy | none | 226,216,800 in the certificate | - | the heavy 3-D storage row blows the like-sized slack up (Gram blocks [10944 4212 x6]); the run hit the 20-minute limit building it, the 'Q' control did not run |
