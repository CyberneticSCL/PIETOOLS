# Shivakumar, "SDP Size Reduction in Polynomial Optimization Problems" (11/17/2026): what it gives for degree sizing (read 11/18/2026)

Source: C:/Users/mpeet/Downloads/newtonPolytope (2).pdf, 22 pages (Sec. 1 main results; App. A-C older attempts).

## Orientation against the toolbox
The document writes Q x(s) = int K_Q(s,theta) x(theta) dtheta: theta is the input variable, s the
output variable (integrated out in Q*Q). The toolbox's lift (sopvar_lift_notes, Tab. slots):
the lift degree D is the degree in the basis operator's own (input) variable, the weight degree
w the half-degree of the weight in theta (the output variable). So the document's
"input-variable bound" is a bound on D, its unbounded "output-variable degree" is w. Its
K_Q = A + B 1_{theta>=s} (global part A, upper-triangular part B) is the toolbox's complete
moments (len_com) plus the incomplete moments (len_inc), up to the lower/upper orientation.

## Results, with their context
1. Lemma 1.2 (triangular degree bound). Context: 1-D, [0,1], P self-adjoint and purely
   integral (no multiplier, no R^q block), kernel of degree <= p in each variable, P = Q*Q with
   Q purely integral. Then the triangular part B has total degree d_{s+theta}(B) <= p - 1.
2. Theorem 1.3 and 1.4 (input-variable bound, reduced basis). Same context. If any
   factorization P = Q*Q exists, one exists whose global part has input degree <= p (and output
   degree <= the output degree c of the given factorization, not bounded by p), and whose
   triangular part has total degree <= p - 1; a Gram form over that reduced basis represents P.
3. Theorem B.4 (Newton polytope, the classical polynomial result) with Lemma B.1 (Sylvester):
   a Gram basis monomial outside half the Newton polytope of the target polynomial has a zero
   diagonal entry, hence a zero row, and can be dropped.
4. Negative results (App. B.3): neither leading-order degree function tried for 1-D PI
   operators (a+b+1, or b) satisfies the non-cancellation property; counterexamples given. So no
   Newton-polytope rule for the integral bases follows, and no bound on the output-variable
   (weight) degree: the document leaves c free.
5. App. C (facial reduction): one paragraph, no result.

## What transfers to the toolbox
- The lift degree D of a certificate's global (complete-moment) block is bounded above by
  p = Dmax(target) for pure integral targets: with the toolbox's necessary D >= Dmin(target)
  (lift notes, Lemma shape), D lies in [Dmin, Dmax], and a ladder on D is finite. Today the rule
  sets D = Dmin + dD and raises w; the document says w is the unbounded direction, D is not.
- The incomplete block of such a certificate needs lift + weight total degree <= p - 1: a
  target-based per-direction cap (the toolbox's cap knob implements theta + s <= cap, but is
  set from Jmax(T) + cap_delta and used in 3-D only). Light's 1-D row (len_inc 2, weight
  half-degree 1) reaches total degree 2, so it certifies targets with p >= 3 at no waste and
  over-parametrizes p = 2 targets by the degree-2 monomials of the incomplete block.
- Newton-polytope pruning of the multiplier (evaluation) block is sound in the separated form:
  the delta part of A*MA comes only from multiplier x multiplier products (composition table,
  M_a* M_b -> M_{a+b}), so the multiplier part of the target is a plain polynomial SOS in the
  basis (weight monomials x evaluation monomials), Theorem B.4 applies to it, and Sylvester
  removes the whole basis element. Gain only when the target's multiplier has a sparse
  support (missing degrees, parity). poscopvar_direct's options.include can carry the mask.
- Per-monomial pruning of the integral bases is NOT justified (App. B.3), consistent with the
  10/2026 finding that possopvar degree pruning is invalid.

## What does not transfer, or needs a proof
- Targets with a multiplier or R^q block: the cross terms (evaluation x integral give kernels
  S(s) Z(theta) on a triangle) can cancel the leading homogeneous term of B*B, which the proof
  of Lemma 1.2 excludes by assuming Q purely integral. Most LPI targets have multipliers.
- Two or more variables, vector intervals: not treated.
- The free operator Q of the direct stability form (deg_Q): not a certificate, not treated; the
  containment rule deg Q >= deg P + deg T of 11/18 is the same kind of statement.

## Measurements worth running (1-D, workstation)
1. Set the certificate's lift degree to Dmax(target) (no dD, no D ladder) with the weight ladder
   only, on the 1-D library rows: cost and verdicts against the current rule.
2. The triangular cap lift + weight <= Dmax - 1 on the incomplete block of the slack, on pure
   integral targets first (a Poincare-type LPI), then on stability targets with multipliers to
   see whether it loses certificates.
3. Newton pruning of the multiplier block on targets with sparse multipliers (even
   coefficients, DEMO3's Poincare target): basis size and verdict.

## Correction (maintainer, 11/18 02:40)
"The targets that matter for complexity come from the derivative, and these almost never carry
multipliers." Right: the negativity target of the direct form, A'Q + Q'A with Q shaped like T,
and of the storage form, T'PA + A'PT, have no multiplier on the L2 part whenever T is purely
integral (every PDE state; the ODE part of a DDE or coupled PIE enters through R^q blocks, not
a multiplier), and the positivity target of the direct form, T'Q - eppos T'T, is purely integral
too. So Lemma 1.2 and Theorem 1.4 apply directly to the slacks that set m and the Gram sizes;
the multiplier caveat above concerns the storage P and the R^q couplings only. Consequences to
measure first, on the 1-D stability rows: the negativity slack at lift degree D = Dmax (no D
ladder) and the incomplete block capped at lift + weight <= Dmax - 1; then the same for the
positivity slack of the direct form.
