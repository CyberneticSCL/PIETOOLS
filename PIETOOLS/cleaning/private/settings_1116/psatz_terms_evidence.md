# Positivstellensatz terms in the N-D certificates: theory and the measured record

Compiled 2026-11-15 (session *New Posopvar and Duality*) at the request of the maintainer, who asked for the evidence that the single-face terms (s_k - a_k), (b_k - s_k) and the cross-direction products (a face or box factor of s_1 times one of s_2) are necessary, against the theoretical results that indicate they are not, and against the 4^N growth of the full term set. Sources: the GPT handoff *Claude_PIETOOLS_multivariate_2026-11-15.zip*, `sopvar.pdf` Sec. 8 (S. Shivakumar, 2026-11-15), the memory notes and scratch logs of this session, the reports under `testing_codes/private`, and the replies of the sessions *Posopvar structure*, *cuADMM for PIE optimization evaluation*, *Low rank PIE/BM solver assessment* and *Repo cleanup and testing* (Sec. 5). Every number below is quoted from a record named in the row; "inferred" marks a statement that no run supports directly.

Notation. Box Omega = prod_k [a_k, b_k]; faces f_{k,1} = (s_k - a_k)/L_k, f_{k,2} = (b_k - s_k)/L_k; box factor g_k = (s_k - a_k)(b_k - s_k); sigma denotes a (matrix) sum of squares. QM(h_1, ..., h_r) = {sigma_0 + sum_i h_i sigma_i} is the quadratic module of the h_i; the preordering adds every product of distinct h_i. The lift degree D (`mult`) and the weight degree w (`int`) are those of `sopvar_lift_notes` Sec. 3; a term "at w - 1" has its Gram one weight degree below the plain term.

## 1. The two documents

### 1.1 The GPT handoff (2026-11-15)

Claims. (a) P(Q) = A_D^* M_W A_D with W = sigma_0 + g_x sigma_1 + g_y sigma_2 + g_x g_y sigma_3, each sigma_i a full matrix Gram form, is sufficient for positivity. (b) By the 2026 preprint of Dritschel, Klep, McCullough and Volcic (arXiv:2608.23073; abstract: "Factorization for positive semidefinite matrix-valued polynomials over a nonsingular compact affine real surface is established. Corollaries include Fejer-Riesz factorization for bivariate matrix polynomials that take positive semidefinite values") every polynomial matrix positive semidefinite on a rectangle, singular values allowed, has the form (a) at some finite degree; the handoff derives this as a corollary (pull-back to the torus, four-sign averaging). (c) By Scheiderer (Spectrahedral shadows, 2018, cited as Cor. 4.24; not verified here) the cone of bivariate polynomials of degree at most six nonnegative on a rectangle is not a spectrahedral shadow, so no fixed-degree affine PSD description of the positive operators exists. (d) Exact obstruction: on [0,1]^2, g_x g_y is not in QM(g_x, g_y) at any degree (the SOS factors are forced to vanish on the faces and the monomial xy cannot be produced). (e) MOSEK passed 17 of 17 SDP checks at D = (0,0), r = (1,1) with supplied weights; the ablation without the g_x g_y term is infeasible under MOSEK and SeDuMi.

What it does not address. The handoff never uses a single face. Its localizer set is exactly the 2^N product set of the box factors, and every test supplies the weight. It makes no statement about tightness at a fixed degree, and its completeness result (b) carries no degree bound ("not a degree bound from the paper"). Its corollary (b) is the handoff's derivation, not a statement of the cited abstract.

### 1.2 `sopvar.pdf` Sec. 8 (2026-11-15)

The section states the conjecture P = A_{d1,d2}^* M_S A_{d1,d2} and the 1-D knobs, then in 2-D: the box is Archimedean, "so Putinar's representation applies" with S = sigma_0 + g_{x,ab} sigma_1 + g_{y,ab} sigma_2; "experience has shown that this representation does not seem to be enough and that some cases require" the seven-term form with the four single faces; "there is some evidence that we also require the products such as sqrt(g_{x,a}) sqrt(g_{y,b}) - yielding 16 terms in 2D and 64 terms in 3D. Perhaps the strongest theoretical evidence against inclusion of such terms is Arxiv: 2608.23073v1."

Traced below: the first "experience" sentence rests on E1, E5 and E6 (Sec. 4), which compare the faces with the box product and with the 2^N product set at one fixed degree; the "some evidence" for the 16-term set rests on E7 and E8, of which E8 is confounded and E7 is a matched-degree effect explained by Lemma T4. No run has used the 16-term set (Sec. 3).

## 2. What the theory says, with the degree made explicit

T1 (Markov-Lukacs, one variable; cited, to verify: Karlin and Studden 1966, matrix case Dette and Studden 2002). Given w >= 1 and a polynomial matrix W of degree 2w with W(s) >= 0 for all s in [a,b], there exist matrix sums of squares S_0 of degree 2w and S_1 of degree 2w - 2 such that W = S_0 + g S_1 with g = (s-a)(b-s). For degree 2w + 1 the form is (s-a) S_1 + (b-s) S_2, both of degree 2w. Gloss: in one variable the box term at w - 1 loses nothing at the matching degree, and the faces are the odd-degree form.

T2 (Putinar 1993; Schmudgen 1991; matrix case Scherer and Hol 2006; cited, to verify). Given a polynomial (matrix) W with W(s) > 0 (positive definite) for all s in Omega, there exist sums of squares such that W is in QM(g_1, ..., g_N); also W is in QM(f_{1,1}, f_{1,2}, ..., f_{N,2}), since each g_k = L_k^2 (f_{k,1} f_{k,2}^2 + f_{k,2} f_{k,1}^2) (from f_{k,1} + f_{k,2} = 1) is in the quadratic module of the faces. No bound on the degrees of the sums of squares is part of either statement as used in the toolbox. Gloss: as the degree grows without bound, the N + 1 terms {1, g_k} suffice for strictly positive weights; faces and products are not needed in the limit.

T3 (handoff corollary of arXiv:2608.23073, derived, to verify). In two variables, every polynomial matrix positive semidefinite on a rectangle, singular values allowed, lies in the preordering of {g_x, g_y}: the four terms {1, g_x, g_y, g_x g_y}, at some finite degree. Gloss: for weights with zeros the product g_x g_y is part of the complete form; T2 does not cover such weights.

T4 (elementary; proven here, the argument of the handoff's Sec. 6). On [0,1]^2 let W = xy. Then (i) W is in no QM(x, 1-x, y, 1-y) and in no QM(g_x, g_y), whatever the degrees of the sums of squares; (ii) W = (xy)^2 + y^2 g_x + x^2 g_y + g_x g_y, a preordering certificate with multipliers of degree 2; (iii) W = f_{1,1} f_{2,1} with the constant multiplier 1.
Proof of (i). Suppose xy = sigma_0 + x sigma_1 + (1-x) sigma_2 + y sigma_3 + (1-y) sigma_4 as polynomials with sigma_i sums of squares. On the segment x = 0, 0 < y < 1 the left side vanishes and every term on the right is nonnegative, so sigma_0, sigma_2, sigma_3, sigma_4 vanish there; a sum of squares vanishing on x = 0 is divisible by x^2. On y = 0, 0 < x < 1 the same argument makes sigma_0, sigma_1, sigma_2, sigma_4 divisible by y^2. Hence every term on the right has y-degree at least 2 or x-degree at least 2, and the coefficient of xy on the right is 0, while it is 1 on the left. For QM(g_x, g_y) the same restriction gives x^2 | sigma_0, sigma_2 and y^2 | sigma_0, sigma_1, and again no term produces xy. (ii) and (iii) are identities (expand).
Gloss: a weight that vanishes on two faces meeting at a corner is reachable by no quadratic module, at any degree; it needs a product term. The corner product of two faces reaches it with constant multipliers, the box products with multipliers of degree 2. This is exactly the structure of the matched-degree test E7.

T5 (Scheiderer 2018, as cited by the handoff; to verify). In two variables no finite term set at a fixed degree describes the nonnegative polynomials of degree six on a rectangle. Gloss: at a fixed degree every set in Sec. 3 is strictly conservative on some weight, and adding a term never hurts (the cones are nested), so fixed-degree comparisons are empirical.

T6 (interleaving; elementary, proven here; the first inclusion was observed by the cuADMM session). Write C_faces(d) for the cone {sigma_0 + sum_k (f_{k,1} sigma_{k,1} + f_{k,2} sigma_{k,2})} with every sigma a Gram form of degree d (monomials up to degree d), C_box(d, e) for {sigma_0 + sum_k g_k sigma_k} with sigma_0 at degree d and the sigma_k at degree e, and C_prod(d; e_x, e_y; e_xy) for the 2-D product set with the plain term at d, g_x and g_y at e_x, e_y and g_x g_y at e_xy. Then for every d >= 0:
(i) C_faces(d) is contained in C_box(d + 1, d), since f_{k,1} + f_{k,2} = 1 gives f_{k,1} sigma = f_{k,1}^2 sigma + g_k sigma / L_k^2, and likewise for f_{k,2};
(ii) C_box(d, d) is contained in C_faces(d + 1), since g_k sigma = L_k^2 (f_{k,1} f_{k,2}^2 + f_{k,2} f_{k,1}^2) sigma;
(iii) in 2-D the corner term f_{1,i} f_{2,j} sigma lies in C_prod(d + 2; d + 1, d + 1; d), by T4(ii) applied to the corner at hand.
Gloss: the faces and the box factors are each within one Gram degree of the other, in both directions, and a corner product is within two degrees of the product set. No fixed-degree comparison between these families can show more than a one- or two-degree saving; what it shows is the size of the certificate at the degree where it first exists.

Consequence. "Not required" in T2 and T3 means not required as the degree grows without bound. The toolbox runs at w in {1, 2, 3}. Only T1 (one variable) decides the question at the degrees in use. In N >= 2 the theory fixes the asymptotically complete sets ({1, g_k}: N + 1 terms for strictly positive weights; the 2^N products of the g_k for singular weights in 2-D), bounds the degree saved by the faces and corners (T6), and says nothing about which fixed-degree set is tighter or smaller.

## 3. The term sets, their counts, and whether each was ever run

| set | terms (N = 2 / 3) | `poscopvar_direct` form | asymptotically complete for | measured |
|---|---|---|---|---|
| none | 1 / 1 | code 0 | nothing (T1, E11) | yes |
| box product g_1 ... g_N (stock 2-D `psatz = 1`) | 2 / 2 | codes [0 1] | nothing in N >= 2 (it vanishes on every face) | yes |
| Putinar box set {1, g_k} | N + 1 = 3 / 4 | codes [0, 2N+2+d] | strictly positive weights (T2) | once (E1 row {g1,g2}; E7 row AQ) |
| product set {prod_{k in S} g_k} (handoff; "tensor set" of 10/08) | 2^N = 4 / 8 | codes [0, 2N+2+d, 1] in 2-D | all PSD weights in 2-D (T3) | three times (E1 at P = I, E6 twice) |
| faces {1, f_{k,1}, f_{k,2}} ('auto' in N-D) | 2N + 1 = 5 / 7 | codes [0, 3 : 2N+2] | strictly positive weights (T2) | throughout |
| corners ('corners', opt-in since 8659c287) | 2N + 1 + 4 C(N,2) = 9 / 19 | faces plus rows with two entries in {1,2} | as faces | E7, E8, E9 |
| full preordering of the faces | 4^N = 16 / 64 | rows in {0,1,2,3}^N | strictly positive weights (Schmudgen) | never |

The 4^N set that the maintainer's section counts has not been built or solved by any session. The largest sets run are nine terms (E1 "3x3 grid", E8/E9 'corners'). The named sets of `sopvar_settings_proposal.tex` Sec. 4.3 are none, box, faces, corners and full with 3^N - 1 nonzero rows over {0,1,2}^N; the 4^N count adds the box factor as a fourth per-direction choice.

## 4. The measured record

Columns: what was compared, the result, the cost, the record. "Reach" is the largest certified fraction of the exact margin. All solves MOSEK unless stated.

### E1. Generator sets with the Lyapunov operator pinned to the identity (cuADMM session, 2026-09-23)

Plant x_t = Lap x + lam x on [0,1]^2, Dirichlet on four edges, lam* = 2 pi^2; P = I, so the slack certificate alone decides; stock `poslpivar_2d` path (retired), heavy settings_2d, eppos 1e-2, eq_deg = LF_deg + 1 ("Dup 1"). Every term, plain or weighted, uses the same Gram basis (block size 264), so a generator of degree k raises the total degree of its term by k above the plain term: faces +1, box factors +2, the product +4. A certificate requires numerr 0, feasratio > 0.5, rel_b < 1e-6 and a non-trivial point; PSD was not checked independently; the sweep stops at the first failure, so a set that never failed is censored at the top of its grid (0.9999). Record: `testing_codes/private/PIETOOLS_demos/cuadmm/results/2026-09-24/logs/{psatz2d,psatzgrid,psatzdeg,perdir,matchdeg}.log`; scripts under `cleaning/private/definitely/PIETOOLS/PIETOOLS_demos/cuadmm/experiments/`; memory `psatz-generators-2d`; reply of the cuADMM session (Sec. 5).

| set | terms beyond the plain one | m | blocks | reach (first failure) |
|---|---|---|---|---|
| none | 0 | 2880 | 1 | fails at 0.10 (rel_b 1.2e-2, numerr 2) |
| g1 g2 (stock `psatz = 1`) | 1 | 5848 | 2 | fails at 0.10, threshold only (rel_b 1.4e-6, numerr 0) |
| ball (stock `psatz = 2`) | 1 | 4676 | 2 | 0.90 (fails at 0.99, numerr 2) |
| g1 alone | 1 | 4228 | 2 | 0.50 (fails at 0.75, numerr 2) |
| g2 alone | 1 | 4228 | 2 | 0.50 (fails at 0.75, numerr 2) |
| {g1, g2} (Putinar set) | 2 | 4676 | 3 | 0.99; fails once at 0.999 (rel_b 5.8e-5, numerr 2); 0.9999 not run |
| {g1, g2, g1 g2} (the 2^N product set) | 3 | 5992 | 4 | >= 0.99, top of its grid; never run above 0.99 |
| four faces ('linear4') | 4 | 3728 | 5 | >= 0.9999 (censored) |
| four corner products of faces | 4 | 4328 | 5 | >= 0.9999 (censored) |
| 9 products {s, 1-s, s(1-s)}_{s1} x {s, 1-s, s(1-s)}_{s2}, no single-direction term | 9 | 5992 | 10 | >= 0.9999 (censored) |

Degree controls: the product at eq_deg - 1 and eq_deg - 2 fails at 0.10 (numerr 2); the faces at eq_deg - 1 or eq_deg - 2, per direction or uniform, fail at 0.90 with the trivial point, and reach >= 0.9999 only at the full eq_deg. Reading: at this Gram degree the product alone certifies nothing; the Putinar set {1, g1, g2} reaches 0.99 and its only measured separation from the faces and corners is one MOSEK failure at 0.999; the faces have the smallest m of any certifying set (3728 against 4676 for the Putinar set). The earlier memory entry "3 x 3 grid >= 0.50" was a misreading; the log gives >= 0.9999. This is the only head-to-head of {1, g_k} against the faces in the record, and the only run of the 2^N product set at P = I.

### E2. Generality of E1 (cuADMM session, 2026-09-24)

Same harness, only none, the product and the faces run (no Putinar set, no corners), grid starting at 0.25 unless stated. Faces: anisotropy >= 0.99; domain offset 0.90 (fails at 0.99 on status only, rel_b 1.8e-7, feasratio 0.40); mixed Dirichlet/Neumann in s2 >= 0.99 (none reaches 0.10 there); advection b1 = 3: 0.50, the last rung below the analytic cap 0.898 of pinned P, failing at the first rung above it (0.90), so the faces certified every rung the pinned form allows (Repo-cleanup session, correcting the memory line "56% of attainable"); [0,1] x [0,2] raw generators fail at 0.25 (threshold only, rel_b 1.5e-6), normalised >= 0.99 (the repo ships normalised faces); [0,1] x [0,4] fails at 0.10 under both forms. The product fails at 0.25 on the baseline, anisotropic, offset and advection plants (rel_b 3.4e-3 to 0.44, numerr 2) and returns the trivial point on the 1 x 2 domain; none fails at 0.25 on all of them. ODE-coupled plant (first with a nonzero multiplier cell): lam_max none 1, product 2, faces >= 19. Locally positive coefficients (lam grid 1, 2, 5, 10, 15, 19): L0 constant: none fails at 1, product 1 (fails at 2 on status only), faces >= 19; L1 reaction lam(1 - 2 s1) x: none and product fail at 1, faces >= 19; L2 diffusion 0.05 + s1(1 - s1): none and product trivial at 1, faces 10; L3 interior patch: none 10, product 15, faces >= 19. Record: `.../2026-09-24/logs/{genplant,normgen,localcoef,odegen}.log`; memory `psatz-generators-2d` (whose line "none = product = 0.00 on every plant" the cuADMM session now corrects: none reaches 0.10 on mixed D/N and lam 10 on L3).

### E3. Searched Lyapunov operator, light storage (cuADMM session, 2026-09-24/25)

Same plant, `lpisettings('light')`, `eq_use_psatz = [3;4;5;6]` at the full `eq_deg`:

| cfg | lam/lam* | m | blocks | nnz(At) | rel_b | numerr |
|---|---|---|---|---|---|---|
| psatz off | 0.10 / 0.50 / 0.90 | 3456 | [8 424] | 557,492 | 1.5e-3 / 2.08 / 1.00 (X = 0) | 2 / 2 / 0 |
| faces | 0.10 / 0.50 / 0.90 | 4492 | [8 424 x5] | 3,658,932 | 4.2e-8 / 8.4e-9 / 1.1e-8 | 0 |

Psatz off under SeDuMi on the same programs returns PSD points with rel_b 7.3e-6, 1.9e-5 and 1.3e-4 at 0.10, 0.25 and 0.50 (numerr 1): under the 1e-6 rule none certifies, so the memory line "psatz-off reach about 0.25 under SeDuMi" applied a looser reading (cuADMM session, Sec. 5). The faces are the only configuration that certifies at any fraction under the rule, to at least 0.90 (not run above), at +30% m and 6.6x nonzeros. At the heavy tier with P searched (`transfer2d`, eppos 1e-2): psatz off fails at 0.25 with the trivial point (m 3096); the faces reach 0.99 and fail at 0.999 on the threshold only (rel_b 2.7e-6, numerr 0; m 3912). The "2-D stalls at 0.25 lam*" figure of the earlier notes comes from stock heavy with psatz off (`search2d`: 0.25 at m 8950, fails at 0.50), not from the product. Corners and the Putinar set were never run with P searched. Record: `.../2026-09-24/logs/{psz_blocks,psz_blocks2,psz_reach,transfer2d,search2d}.log`; memory `psatz-generators-2d`, `psatz-2d-linear4-verified`; `testing_codes/private/PIETOOLS_demos/cuadmm/results/2026-09-26/REGIME_REPORT.md` Sec. 5 (m = 4492, MOSEK F in 64 s at 0.5 lam*; cuADMM not F after 10k iterations). The container baseline rows `stab2_rd` and `stab2_rd_psz` of the cuADMM suite assemble the same SDP (the face codes are not applied on the container settings) and are evidence of nothing.

### E4. Faces at a point where psatz off already certifies (Low-rank session, 2026-09-25)

Same plant, 0.10 lam*, Dup 3, solver 'best': off m 3456, rel_op 9.46e-7 (SeDuMi, 295 s); faces m 4492, rel_op 5.40e-8 (MOSEK, 1198 s). No change of verdict; residual 17.5x (operator) and 179x (row) tighter; cost 4.1x time, 6.6x nnz. The 0.50/0.90 sweep of that session was stopped and never completed; no low-rank rank per generator set exists. Record: memory `psatz-2d-linear4-verified`; transcript 1a573830 line 4493.

### E5. The 2-D heat benchmark and the Poincare diagnostic: the product alone is certified infeasible, the faces reach the ceiling (Posopvar session, 2026-09-27/29)

All rows MOSEK 11 at tolerance 1e-10; a bracket [lo, hi) has lo certified feasible and hi certified infeasible by a verified Farkas ray; "uncertain" is MOSEK UNKNOWN and is evidence of nothing. 'product' is g1 g2 (the stock `psatz = 1`), never the Putinar set {g1, g2}, which did not exist as codes before 10/08. Record: `PIETOOLS_demos/sopvar_demos/README.md` Secs. 3, 6, 7, 9, 11.3, 11.4; scripts `heatNd_poincare_validate`, `heatNd_validate`, `fq_t2d`; reply of the Posopvar session (Sec. 5).

(a) 2-D Poincare diagnostic, L2[s1,s2], DD x DN, d = 1, Gram block 64 (the d = 1 basis ceiling is 2.0 for the full L2 form, 1.0 per direction; the energy form has the exact value 12.337006):
| form | none | product g1 g2 | four faces |
|---|---|---|---|
| L2 full | [0, 0.032), m 512 | [0, 0.090), m 1152: certified infeasible at 4.5% of the ceiling | 2.0000000 at the ceiling, m 768, 5 x 64 |
| L2 per direction | [0, 0.0072) / [0, 0.0096) | [0, 0.039) / [0, 0.018) | 1.0000000 / 1.0000000 |
| energy full | - | [0.337, 1.349), m 1376: certified infeasible at 10.9% of exact | lo 12.035914 (97.6%), m 920; d = 2: 12.3368080 |
The stock `poslpivar_2d` route gives the same numbers (product 0.00241 / 0.00132 / 0.00048; faces 1.99994 / 0.99997 / 0.99997). These are the only Farkas-certified failures of a term set in the whole record, and both concern the product alone.

(b) Face subsets on the same L2 target at d = 1 (`fq_t2d`, 09/28; the stock face codes and the container codes give identical SDPs and identical lo): all four faces lo 1.9999443; the upper face of s1 with the lower face of s2 lo 0.0240957; one face alone lo 0 (no positive lam certified in 30 bisection attempts; not a certified infeasibility). So at d = 1 every face is needed on this target, in the sense of what certifies.

(c) Heat LPI, Cor. 35, u_t = u_{s1 s1} + u_{s2 s2} + kappa u, DD x DN, kappa* = 12.337006. No psatz: certified infeasible at 0.048 kappa* (d = 0) and 0.386 kappa* (d = 1). Faces at bench (d = 0): [12.183790, 12.183978), m 1242, K.s 48 x 5 + 76 x 5, nx 40,409; heavy: lo 12.336941, m 4282. Ablations of bench d = 0: the product in place of the faces is certified infeasible at 0.5 kappa* while the faces certify 0.988 kappa*; dropping the Neumann face or the s2 Dirichlet face gives MOSEK UNKNOWN (rel_b 3.5e-3 to 7e-3, no ray) at kappa = 0, 0.1 and 0.5 kappa*, which is not a proof that the face is needed; faces on Q only is certified infeasible. Table 1 of the paper does not reproduce from the printed listing without the faces (maintainer, 09/27: a version difference). Also: `heatNd_lpi.m` lines 72-99, `copquadvar.m` header MMP 09/27/2026, memory `heatNd-benchmark`.

(d) 3-D heat at bench with the faces (7 terms per operator, Gram 288 and 552, m 31,770, 2.71M variables): kappa in [14.60, 14.70) of 14.804407, 499-921 s, 20-21 GB. No term ablation in 3-D. 3-D Poincare per-direction targets with only the two faces of the tested direction: certified at 0.99 of the prediction (L2, m 20,480, 3 x 512, about 600 s, 9.7 GB; energy at 0.99 and 0.999); the sum of the three per-direction certificates certifies the all-faces L2 problem at 2.97 without a solve; the 1-D certificate tensored with the other directions needs only the g of its own direction (derived).

### E6. The 2^N product set against the faces at matched per-direction degree (this session, 2026-10-08)

The "tensor set" {plain at w, g_d at w - 1 in direction d alone, g_1 g_2 at w - 1} is the handoff's four-term form; it was added as codes 2N+2+d and compared with the faces at w - 1 (same per-direction degree box as the plain term) and at w.

| problem | faces at w | faces at w - 1 | product set at w - 1 | product set at w | product pair |
|---|---|---|---|---|---|
| 2-D coercive KYP io2, (2,2) | gamma 0.171718, nx 557,845, 57 s | 0.171722, nx 200,425, 16 s | 0.171720, nx 233,280, 17.6 s (cap 6: 0.171775, 197,720) | - | (2,2), (3,2), (2,3): UNKNOWN; (3,3): 0.171777 at nx 1,427,142, 491 s |
| 2-D heat DD x DN, R (1,1), Q (2,2), d = 0 | nx about 1.04M (over budget) | kappa 12.336731, nx 339,209, 16 solves | 12.326178, nx 395,785, 15 solves, 2x time | nx 762,889, two solves in budget | nothing |

Reading: on the KYP slack the product set equals the faces at w - 1 within 2e-6 in gamma at 1.16x the variables; on the heat benchmark it certifies a gap of 8.8e-4 against 2.2e-5 for the faces at 1.17x the variables and twice the time. The product pair (plain + g1 g2) needs one more degree on KYP and certifies nothing on heat. This made faces at w - 1 the N-D default (commit 96ef81f6). Record: memory `certificate-dual-structure`, `hinf-slack-sizing`, `heatNd-degree-tailoring`; `sopvar_lift_notes.tex` Sec. 10, Table 3.

In one variable, same date: faces at w against the product pair at w - 1 (the Markov-Lukacs form) give the same gamma at 2.2x the variables (io1: nx 2019 vs 898); at equal w the faces are tighter (DD heat Q(2,1): 9.8690 vs 9.8048; (1,1) KYP 0.186930 vs 0.193093) because the face cone at w contains the box cone at w - 1 (T2). Record: memory `hinf-slack-sizing`, `heatNd-degree-tailoring`.

### E7. Corner products at the matched degree (this session, 2026-11-12)

<x, P x> = int int s1 s2 (x1(s1) + x2(s1,s2))^2 on L2[s1] x L2[s1,s2], [0,1]^2: the weight is the corner product s1 s2 of T4. Certify P + eps I in the cone of each set at Gram degree D and minimise eps (eps* = 0 means exact). N = total Gram size. Record: scratch `psatz_mix/pm_run.m`, `pm_run2.m` and logs; memory `psatz-corner-products`.

| set | D = 0 | D = 1 | D = 2 |
|---|---|---|---|
| plain only | infeasible | 3.6e3 (numerr 2) | 18.9 (numerr 2) |
| plain + box product | infeasible | 9.4e2 (numerr 2) | - |
| plain + {g1, g2} (Putinar set) | infeasible | 8.9e-3 | - |
| plain + four faces | infeasible | 8.9e-3 (N 600) | 2.8e-4 |
| plain + faces + {g1, g2} | infeasible | 8.9e-3 | - |
| plain + faces + box product | infeasible | 2.4e-8 | - |
| plain + faces + {g1, g2} + box product | infeasible | 2.0e-8 (N 960, 58,080 dvars) | 1.4e-8 |
| plain + corner product s1 s2 alone | 3.1e-10 (N 24, 156 dvars) | 5.5e-9 | 5.5e-10 |
| plain + four corner products of faces | 1.9e-8 (N 60) | 1.3e-8 | - |
| all codes + eight mixed rows | 1.3e-8 (N 192) | 5.4e-9 | 3.4e-9 |

Reading, by T4: the quadratic modules (faces; {g1, g2}; both) never reach the weight and converge slowly (8.9e-3 at D = 1, 2.8e-4 at D = 2); the sets containing a product term are exact one degree later (D = 1, Gram 960) and the corner product is exact at D = 0 (Gram 24). The corner products save one degree on a weight that vanishes on two faces meeting at a corner; they do not do anything the product set cannot do at D + 1.

### E8. Corner products on the edge-coupled 1-D/2-D heat cascade S3, and the storage confound (this session and the Repo-cleanup and mixed-domain sessions, 2026-11-12/13)

`PIETOOLS_stability_sop` form 'P', mixed-domain system S3 (1-D edge state coupled to a 2-D heat state, reaction r, r* = pi^2), six-step bisections. Record: scratch `psatz_mix/exec/reach_S3_{A,B}_light.log`, `t_exec*.log`; `testing_codes/private/converters/tests/mixed_domain/MIXED_DOMAIN_REPORT.md` "Beyond light"; memory `psatz-corner-products`.

| storage | slack terms | reach | variables, time per solve |
|---|---|---|---|
| light stock (settings_2d) | faces | 0.4844 r* (0.5 unknown at dw 0 and 1) | 318k / 508k, 13-48 s |
| light stock | faces + 4 corners | 0.6719 r* (0.5 at dw 0; 0.6875 unknown) | 572k / 762k, 38-100 s |
| heavy stock | faces | 0.984 r* (= bisection resolution) | 1.02M, 95 s |
| heavy stock | faces + corners | 0.984 r* (no gain) | 1.59M |
| light, N-D default storage (DEG omitted, 'def' routing) | faces | 0.9 r* certified at dw 0 (eta 1.2e-10; S8 eta 5.6e-11 at 335,104); at 1.1 r* infeasible | 318,004, 14-22 s |

The N-D-default rows were measured first by the Repo-cleanup session on 11/12 (019b383e export, light, MOSEK; records `scratchpad 874c85cb/p3/{probe,sweep}.out`) and again by this session's routing runs (`psatz_mix/route/rt_route.log`); that session infers, from the equal slack sizes [252 x 5], that the settings_2d storage sets the light-tier ceiling, and marks the inference as such.

The same holds for S8 (two edge states): light stock faces reach in [0.575, 0.6]; corner rows certify 0.5 at dw 0 but not 0.9; heavy faces certify 0.9; the N-D default storage with faces certifies 0.9 at 335k. S17 is certified by no set at any tier. So the gain of the corner rows at light was a property of the 2-D stock storage of the Lyapunov operator, not of the slack: with a different storage the faces alone reach 0.9 r* in fewer variables than the corner rows. The corners were made opt-in on 11/13 (8659c287) on the maintainer's instruction.

### E9. Corner products on single-space and in-domain problems (this session, Sol, 2026-11-12)

2-D heat benchmark DD x DN: faces and faces + corners both certify 12.183764 at the bench preset (d = 0 and d = 1), corners at 1.7x the solve time; at heavy both certify 12.336731, 375 s against 623 s. In-domain and ODE couplings (S6, S7, S9) certify at dw 0 with the faces. 3-D: corners never run (19 terms). Record: memory `psatz-corner-products`; `lpi_programming_sopvar/README.md` lines 228-232.

### E10. The handoff's own runs (2026-11-15)

Eight supplied weights and their double-Volterra operators at D = (0,0), r = (1,1): MOSEK 17 of 17 within residual 6.8e-9; SeDuMi 15 of 17 (two residuals 1.7e-7 and 2.9e-7); the ablation of the g_x g_y term on the weight g_x g_y is infeasible under both solvers, as T4 predicts; the single faces were not run. Record: `results/VALIDATION_SUMMARY.txt`, `results/{mosek,sedumi}/matlab_results.json` of the handoff.

### E11. One variable: the box term is required and the faces are not (several sessions)

Poincare, psatz off: gamma H1'H1 - H2'H2 is in the plain cone at no degree and no gamma (proved 2026-09-28; memory `poincare-psatz-off-infeasible`). P = I, 1-D reaction-diffusion, SeDuMi, light, eppos 1e-2 (`psatz1d.log`, `poinc1d.log` of the cuADMM session): none 0.00 at Dup 1, 2, 3; the box term s(1-s) 0.50 at Dup 1 and 0.9999 at Dup 2; the two faces {s, 1-s} 0.10 at Dup 1 and 0.00 at Dup 2; box term plus faces 0.50 at Dup 1, no gain. Confound noted by the Repo-cleanup session: every psatz term in that log sits one degree below the base (the stock offset), and `psatzdeg.log` shows that the same offset removes the reach of the faces in 2-D, so the 1-D face rows there are not at matched degree; the Posopvar rows below are. Searched P at the light preset: the box term changes the verdict at 0.10 and 0.50 lam* where psatz off returns X = 0 (Low-rank session reply, `pielr_settings.m` header lines 49-58). A face term at weight 0 certifies nothing in the 1-D Dirichlet case (`get_lift_degs` header 10/08). Posopvar session, 1-D Poincare and heat (README Sec. 11.3, 09/27-28, both terms at the full Gram degree): DD L2 at d = 1, 2, 3 the box term gives 1.0000000 / 9.8692933 / 9.8696044 (m 24 / 51 / 88) and the two faces 1.0000000 / 9.8694053 / 9.8696044 (m 20 / 45 / 80); DD energy d = 1: 9.7013012 against 9.7165200; DN: identical. Heat 1-D: the product with a tensor P at d = 1 gives [9.869552, 9.869628), the same as heavy with the faces. Every 1-D tier uses the Markov-Lukacs pair and no face. In one variable the record agrees with T1: at d >= 2 the two forms give the same value; at d = 1 the faces at the full degree are a slightly larger cone (T6 (ii)) with fewer rows, and the Markov-Lukacs pair at (w, w - 1) is the complete form.

## 5. Replies of the other sessions

- *Low rank PIE/BM solver assessment* (received 2026-11-15): the searched-P table (E3) is the cuADMM session's; its own measurement is E4 (one point, no verdict change, residuals 17.5x / 179x tighter); its 0.50 and 0.90 sweep never completed; no low-rank rank per generator set exists; the 1-D box-term table of E11. It notes that `set2d_deg.m` (lines 18-21) states that a single `poslpivar_2d` call with `psatz = 1` builds the weighted term instead of the plain one; `lpi_ineq_2d` adds the two calls (`Deop + Deop2`, header lines 28-38), and the block counts of E1 show the plain term present in every row, so this does not affect E1.
- *cuADMM for PIE optimization evaluation* (received 2026-11-15): the full P = I tables (E1, E2), the searched-P runs (E3) and the 1-D runs (E11), with the record paths given there; the corrections to its memory note (3 x 3 grid >= 0.9999, not >= 0.50; "none = product = 0.00 on every plant" false; the SeDuMi "0.25" a looser reading); its 09/23 "structural gap" text ("the Putinar/Schmudgen certificate on a box needs each generator separately, and the product-only cone is strictly smaller ... a strong candidate for why 2-D stalls at 0.25") labelled by itself as inferred and never proved, with the 0.25 figure traced to psatz-off stock heavy; the observation that a face certificate at Gram degree d lies in the box-factor quadratic module at d + 1 (T6 (i)), so the measured face-versus-box gap is a fixed-degree effect consistent with the theory. It states that nothing in its record proves any generator set infeasible: every failure is a solver status or a residual threshold.
- *Posopvar structure* (received 2026-11-15): the 2-D Poincare diagnostic with Farkas-certified brackets, the face-subset runs, the heat ablations, the 3-D bench and per-direction runs (E5) and the 1-D comparisons (E11). It states that its 2-D data test the box product g1 g2 only, never the quadratic module {g1, g2} nor the preordering, and that it never measured corners or any term set for the H-infinity or H2 problems.
- *Repo cleanup and testing* (received 2026-11-15): re-read every number of E1-E3 from the committed logs and confirms them; one measurement of its own, the N-D-default-storage rows of E8; four corrections to the memory note (the 3 x 3 grid; the omitted preordering row {g1, g2, g1 g2} at 0.99, top of its grid; "from 0 to 0.9999" compares against none and the product only; the advection "56%" line); the stock-offset confound of the 1-D face rows (E11). Its summary: in the records found, {g1, g2} never certifies nothing on a stability, KYP or heat LPI where faces or corners certify; every "certifies nothing" row is g1 g2 alone or no psatz; {g1, g2} was run on an LPI only at one degree (E1), and its only failure anywhere is the constructed kernel of E7 at D <= 1.

## 6. Assessment

1. No measurement supports the 4^N set. It has never been built. The sets that were run and gained something have 2N + 1 terms (faces) or, in one matched-degree test, one corner row; the 'corners' set adds 4 C(N,2) terms, quadratic in N.
2. The theory does not contradict the measurements, because it is asymptotic in the degree (T2, T3) while every measurement is at w in {1, 2, 3}, and T5 says no fixed-degree set is complete in N >= 2. The one fixed-degree theorem, T1, is in one variable, where the record agrees with it (E11, E6 1-D rows).
3. Faces against the box sets, at fixed degree, in 2-D: the product g1 g2 alone certifies nothing on every plant measured (E1, E2, E5, E6), and a plain term plus that product is exactly the stock 2-D `psatz = 1`; the only Farkas-certified failures of a term set in the record are its two (E5 (a) and (c)), and no run has certified the Putinar set or the 2^N product set infeasible anywhere; the Putinar set {1, g1, g2} reaches 0.99 and is separated from the faces (>= 0.9999) by one MOSEK failure at 0.999, at 1.25x the rows (E1, the one head-to-head); the 2^N product set reaches the top of its grid at P = I (E1), equals the faces on the KYP slack and is 40x looser in gap on the heat benchmark, at 1.16-1.17x the variables (E6). By T6 the faces can save at most one Gram degree over the box sets; the record says that at the degrees in use they do, on the heat benchmark (E5, E6) and on the pinned-P plants (E1, E2), and that they are the smallest certifying set there. In 3-D the faces (7 terms) against the product set (8 terms) is unmeasured.
4. Corner products: necessary only in the matched-degree sense of T4 and E7, on weights vanishing on two faces that meet at a corner, where they save one degree over the product set; no gain on the heat benchmark at any tier (E9), none at heavy on S3 (E8); the light-tier gain on S3 and S8 was a storage artefact (E8). Opt-in is consistent with the record.
5. The sentence in `sopvar.pdf` Sec. 8 that Putinar's form "does not seem to be enough" is supported as a fixed-degree tightness statement (E1, E5, E6) and not as an insufficiency statement; the sentence that "some cases require" the 16 products has no measurement behind it. The nearest run, the nine two-direction products at P = I (E1), matched the faces at 1.6x the rows, and the only case where a corner row did what no other set did at the same degree is the constructed weight of E7.
6. Cost is not the term count. The faces have the smallest m in E1 (3728 against 5848 for one product term) because they add degree-1 generators; the right cost measures are Gram entries, rows and nnz (E3: +30% m but 6.6x nnz).

## 7. What would settle the question

(a) A head-to-head of {1, g_k} (N + 1), the 2^N product set, and the faces (2N + 1), each at matched per-direction degree, on the four 2-D problems with substantive gaps (heat DD x DN, KYP io2, S3 and S8 with the N-D default storage, the PDE 48 estimator) and on the 3-D heat benchmark: reach, gap and variables. E1 and E6 are the only such comparisons and cover three problems. (b) Reach against w for each set on one problem, to measure the degree at which {1, g_k} catches the faces (the asymptotic claim of T2); the cones are nested in w, so the curves are monotone. (c) For the corners, a problem whose certificate must vanish at a corner at the margin, if one exists among the executives; none has been identified.

## 8. Addendum: the measurements of 2026-11-15 (autonomous run on the maintainer's word)

The plan of `psatz_faces_test_plan.md` was run on the workstation (results, decisions and the harness in `sopvar/psatz_faces_results/`; read `SUMMARY.md` there first). With the lift sized from the target and the box factors lowered in their own direction only, the Putinar set {1, g_k} gives the faces' verdicts on every pinned-P plant of E1 and E2 at 68-75% of the variables, so H0 of the plan holds and the E1 separation was the degree structure of that run. At w >= 2 in two variables every set (box, 2^N products, faces at either offset rule, corners) ties on the heat benchmark at rtol 1e-4, on the io2 gain, on S3 and S8 and on the library rows; the scalar-rule faces are the smallest certifying configuration there (box +5% to +17%). What survives of the faces' necessity is the odd-degree half step: with the plain block at weight w the faces at offset 0 certify strictly more than any box configuration at the same w (Farkas-separated in 1-D, 2-D and on the pinned plant), certify the PDE 48 estimator one weight degree earlier (0.00451 at 110,135 variables against today's 0.00618 at 195,815), and in three variables at the bench caps are the only set that certifies 0.985 kappa* (the box set needs 29.6M variables at the next cap, not solvable here). Corners changed no verdict in any case (the S3 gain of E8 was the storage); the product g1 g2 alone remains the only Farkas-infeasible set; the 4^N set was never needed. Sections 6 and 7 above stand, with item 3 sharpened: the faces' fixed-degree advantage is the half step, present at the lowest weight and in 3-D, absent at w >= 2 in 2-D.
