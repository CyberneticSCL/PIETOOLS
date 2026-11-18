# Coupled-domain linear PDE candidates for the PIETOOLS mixed-domain suite

Literature search, 2026-11-17, about 45 minutes. No computation and no repository edits.

Provenance labels used below:
- **[read]**: the equations or the value were read in the paper (full text, or the arXiv HTML).
- **[abstract]**: confirmed from the abstract or a catalogue record only.
- **[derived]**: not in the paper. It is a separable reduction worked out here and needs checking before it is used as a reference value.
- **[unverified]**: from memory or a secondary citation, and not confirmed.

Spaces use the vocabulary of the brief: R^q, L2[s1], L2[s2], L2[s1,s2], 3-D. Two rectangles that share an edge are mapped onto one unit square by the reflection s1 -> -s1. The interface then becomes the edge s1 = 0 of both 2-D states. A boundary condition that couples the traces of two 2-D states on the same edge is assumed expressible but has not been checked against the converter.

---

## 1. Wave-heat transmission on two rectangles (Batty, Paunonen, Seifert 2019)

- **Equations [read]:** The wave occupies Omega_- = (-1,0) x (0,1) and the heat equation occupies Omega_+ = (0,1) x (0,1).
  - u_tt = Delta u on Omega_-.
  - w_t = Delta w on Omega_+.
- **Boundary conditions [read]:**
  - u(-1,y) = 0 and u(x,0) = u(x,1) = 0.
  - w(1,y) = 0 and w(x,0) = w(x,1) = 0.
- **Interface at x = 0, for y in (0,1) [read]:** u_t(0,y,t) = w(0,y,t) and u_x(0,y,t) = w_x(0,y,t).
- **Energy:** E = (1/2) int (|grad u|^2 + |u_t|^2 + |w|^2), and E is non-increasing.
- **Known result [read]:** For classical solutions, E(t) = o(t^{-2/3}), and the rate is sharp. The semigroup is therefore not exponentially stable, because trapped rays run parallel to the interface.
- **Spaces:** two 2-D states, the wave as (u, u_t) and the heat as w, on L2[s1,s2], or three 2-D states after the reflection.
- **Coupling mechanism:** a two-way interface between two 2-D states. A trace of each state, and the normal derivative of each, enters the boundary condition of the other on the shared edge.
- **Gaps addressed:** (a) two coupled 2-D states; a 2-D wave; a parabolic-hyperbolic interface.
- **Threshold:** There is no parameter threshold. The result is qualitative. A correct method must certify non-increase of the energy, which is a non-strict Lyapunov inequality, and must fail to certify any exponential decay rate.
- **Feasibility:** Linear, constant coefficients, on boxes. It needs the reflection and a boundary condition that couples two 2-D states on one edge. The wave needs the first-order rewrite (u, v = u_t), so the interface condition reads v(0,y) = w(0,y).
- **Citation:** C. J. K. Batty, L. Paunonen, D. Seifert, "Optimal energy decay for the wave-heat system on a rectangular domain," SIAM J. Math. Anal. 51(2):808-819, 2019. DOI 10.1137/18M1195796; arXiv:1803.05405.

## 2. Wave-heat transmission in 1-D (Batty, Paunonen, Seifert 2016; Zhang, Zuazua 2004)

- **Equations [read]:**
  - u_tt = u_xx on (-1,0).
  - w_t = w_xx on (0,1).
  - At x = -1: u_x(-1,t) = 0. Section 5 of the paper uses u(-1,t) = 0 instead.
  - At x = 1: w(1,t) = 0.
  - At x = 0: u_t(0,t) = w(0,t) and u_x(0,t) = w_x(0,t).
  - dE/dt = -int_0^1 |w_x|^2.
- **Known result [read]:** E(t) = o(t^{-4}) for classical solutions, and the rate is optimal. The resolvent grows as O(|s|^{1/2}) on the imaginary axis.
- **Dirichlet case [read]:** The spectrum is the set of nonzero lambda solving sqrt(lambda) sinh(lambda) cosh(sqrt(lambda)) + cosh(lambda) sinh(sqrt(lambda)) = 0. It lies in the open left half-plane but approaches the imaginary axis.
- **Zhang and Zuazua 2004 [abstract]:** the same type of linearized 1-D fluid-structure model, with a sharp polynomial decay rate from a spectral analysis.
- **Spaces:** L2[s1] states only (wave as (u, u_t), plus heat), after the reflection.
- **Coupling mechanism:** a point interface, with traces in each other's boundary conditions.
- **Gaps addressed:** (a) parabolic-hyperbolic coupling. This is the 1-D precursor of item 1 and is the cheaper first test, in line with test-1-D-before-2-D.
- **Threshold:** There is no parameter threshold. The system is strongly but not exponentially stable, and the eigenvalue equation above gives an exact spectral check.
- **Feasibility:** expressible today, as 1-D states with point coupling.
- **Citations:**
  - C. J. K. Batty, L. Paunonen, D. Seifert, "Optimal energy decay in a one-dimensional coupled wave-heat system," J. Evol. Equ. 16:649-664, 2016. arXiv:1505.04951. The journal DOI is 10.1007/s00028-015-0316-y [unverified].
  - X. Zhang, E. Zuazua, "Polynomial decay and control of a 1-d hyperbolic-parabolic coupled system," J. Differential Equations 204:380-438, 2004. DOI 10.1016/j.jde.2004.02.004. Volume and pages are [unverified].

## 3. Wave with an acoustic boundary condition on one edge (Beale 1976)

- **Equations [unverified, standard form of the model]:**
  - phi_tt = c^2 Delta phi in the domain.
  - On the acoustic edge Gamma_1, each boundary point is a damped spring with normal displacement delta(x,t): m delta_tt + d delta_t + k delta = -rho phi_t, with delta_t = d phi / d nu.
  - On the rest of the boundary: d phi / d nu = 0 or phi = 0.
- **Known result:**
  - Beale 1976 shows that the generator has essential spectrum and that the system is not uniformly (exponentially) stable [abstract and secondary citations].
  - Munoz Rivera and Qin 2003 prove polynomial decay for the energy when an acoustic condition holds on part of the boundary [abstract]. Their domain is a smooth bounded domain, not a rectangle.
- **Spaces:** a 2-D state (the wave, first-order form) plus a 1-D edge state (delta, delta_t) in L2[s1]. The edge state has no derivative in s1.
- **Coupling mechanism:** a two-way edge coupling. The edge state is the normal-derivative datum of the 2-D state, and the 2-D trace phi_t drives the edge equation. This is a Wentzell-type coupling without tangential diffusion.
- **Gaps addressed:** (a) a 2-D wave coupled to a 1-D state; a relative of S4 and S25.
- **Threshold:** None is published. The result is qualitative (no exponential stability), and a box version has not been found.
- **Feasibility:** expressible on the unit square with the acoustic edge s2 = 0, given the first-order wave rewrite. The essential spectrum makes this a stress test: a certificate of exponential stability would be wrong.
- **Citations:**
  - J. T. Beale, "Spectral properties of an acoustic boundary condition," Indiana Univ. Math. J. 25(9):895-917, 1976. DOI 10.1512/iumj.1976.25.25071 [unverified].
  - J. E. Munoz Rivera, Y. Qin, "Polynomial decay for the energy with an acoustic boundary condition," Appl. Math. Lett. 16(2):249-256, 2003. DOI 10.1016/S0893-9659(03)80039-3.

## 4. Linear Turing system: two 2-D states coupled in-domain

- **Equations [standard; unverified for this exact statement]:** On (0,L1) x (0,L2), with homogeneous Neumann conditions on all edges:
  - u_t = D_u Delta u + a11 u + a12 v.
  - v_t = D_v Delta v + a21 u + a22 v.
- **Threshold [derived from the standard mode analysis]:** The cosine mode with wavenumber k^2 = (m pi/L1)^2 + (n pi/L2)^2 is stable iff two conditions hold:
  - tr_k = (a11 + a22) - (D_u + D_v) k^2 < 0;
  - det_k = (a11 - D_u k^2)(a22 - D_v k^2) - a12 a21 > 0.
  - The system is stable iff both hold for every admissible pair (m,n), including (0,0).
  - The critical diffusion ratio is therefore exact once the domain is fixed, as the smallest D_v/D_u for which some admissible k^2 gives det_k <= 0.
- **Spaces:** two 2-D states in L2[s1,s2].
- **Coupling mechanism:** in-domain multiplier coupling between two 2-D states.
- **Gaps addressed:** (a) a 2-D state coupled to another 2-D state in-domain, with an exact threshold.
- **Feasibility:** expressible today, with Neumann conditions only. The (0,0) mode makes the threshold depend on the kinetic matrix alone. A Dirichlet variant uses sine modes instead.
- **Citations:**
  - A. M. Turing, "The chemical basis of morphogenesis," Phil. Trans. R. Soc. Lond. B 237:37-72, 1952. DOI 10.1098/rstb.1952.0012 [unverified DOI].
  - J. D. Murray, *Mathematical Biology II: Spatial Models and Biomedical Applications*, 3rd ed., Springer, 2003, ch. 2 (linear Turing conditions on a finite domain) [unverified chapter].

## 5. Two-layer heat conduction (composite medium), extended to a rectangle pair

- **Equations [abstract for the 1-D multilayer problem]:**
  - u1_t = k1 u1_xx on (0,l1) and u2_t = k2 u2_xx on (l1,l).
  - At the interface: imperfect contact, with flux continuity k1 u1_x = k2 u2_x and a jump condition k1 u1_x = H (u2 - u1).
  - Outer conditions: Dirichlet or Neumann.
  - Hickson et al. give exact eigenfunction-expansion solutions. The decay rates are roots of a transcendental equation [abstract].
- **2-D extension [derived]:** Put each layer on a rectangle, with the same Dirichlet or Neumann conditions in s2. The sine or cosine modes in s2 decouple the problem, so each mode reduces to the 1-D two-layer eigenproblem with an added reaction -k_i (n pi)^2. Adding a reaction term lambda u_i gives a stability threshold lambda* = -(the first eigenvalue of the transmission operator), computable from the transcendental equation.
- **Spaces:** two 2-D states in L2[s1,s2] after reflection, or 1-D states only for the 1-D version.
- **Coupling mechanism:** a two-way interface (flux and Robin jump) between two states on the shared edge.
- **Gaps addressed:** (a) two coupled 2-D states; different coefficients per subdomain, which gives piecewise-constant coefficients without needing polynomial ones.
- **Threshold:** semi-analytic. It is not stated in the paper as a stability threshold [derived].
- **Feasibility:** expressible, given a boundary condition that couples two 2-D states on one edge.
- **Citation:** R. I. Hickson, S. I. Barry, G. N. Mercer, "Critical times in multilayer diffusion. Part 1: Exact solutions," Int. J. Heat Mass Transfer 52:5776-5783, 2009. DOI 10.1016/j.ijheatmasstransfer.2009.08.013. Volume and pages are [unverified].

## 6. Wave with partially delayed boundary control (Xu, Yung, Li 2006)

- **Equations [abstract; exact placement of the gain unverified]:**
  - w_tt = w_xx on (0,1), with w(0,t) = 0.
  - Boundary control at x = 1, split into an undelayed part of weight mu and a delayed part of weight (1 - mu): w_x(1,t) = -k [ mu w_t(1,t) + (1-mu) w_t(1,t-tau) ].
- **Delay as a transport state:** z(t,rho) = w_t(1, t - tau rho) on rho in [0,1], with tau z_t + z_rho = 0 and z(t,0) = w_t(1,t).
- **Threshold [abstract]:**
  - For mu > 1/2, a gain k exists that makes the closed loop exponentially stable for every delay tau.
  - For mu = 1/2, the closed loop is at most asymptotically stable.
  - For mu < 1/2, the closed loop is unstable.
  - The closed loop is a Riesz-spectral system.
- **Spaces:** a 1-D wave in L2[s1] (first-order form) and a 1-D transport state in L2[s2].
- **Coupling mechanism:** a point trace of the wave is the inflow of the transport state, and the outflow enters the wave's Neumann condition.
- **Problem case addressed:** (b2) delay as transport with a sharp analytic threshold (mu* = 1/2). It is the boundary-delay counterpart of S18.
- **Feasibility:** expressible today.
- **Citation:** G. Q. Xu, S. P. Yung, L. K. Li, "Stabilization of wave systems with input delay in the boundary control," ESAIM: COCV 12(4):770-785, 2006. DOI 10.1051/cocv:2006021.

## 7. Wave with delayed boundary or interior damping (Nicaise, Pignotti 2006; Datko, Lagnese, Polis 1986)

- **Equations [read in Fridman and Orlov 2009, Sec. 5]:**
  - z_tt = a z_xx - mu0 z_t(t) - mu1 z_t(t - tau(t)) - a0 z(t) - a1 z(t - tau(t)) on (0,pi), with Dirichlet ends.
  - Nicaise and Pignotti treat the multi-dimensional case with the delayed damping either interior or on part of the boundary [abstract].
- **Threshold [read, Remark 6 of Fridman and Orlov]:**
  - With a = 1 and a0 = a1 = 0, constant delay, and mixed Dirichlet-Neumann conditions, Nicaise and Pignotti show stability for 0 <= mu1 < mu0.
  - If mu1 >= mu0, there is a sequence of arbitrarily small delays that destabilize the system.
  - Fridman and Orlov, Corollary 2: exponential stability holds for all delays with tau_dot <= d < 1 if mu1^2 < (1-d) mu0^2.
- **Datko, Lagnese, Polis 1986 [secondary citations]:** A 1-D boundary velocity feedback that stabilizes the delay-free wave is destabilized by arbitrarily small delays.
- **Spaces:** L2[s1] for the wave, and L2[s1,s2] for the delayed velocity z_t(x, t - tau rho) as a transport in rho. In 2-D (Nicaise and Pignotti), the wave is L2[s1,s2] and the interior delay is 3-D.
- **Coupling mechanism:** delay as transport on an added dimension, interior (item 7a) or at an edge (item 7b).
- **Problem case addressed:** (b2). The threshold mu1 = mu0 is delay-independent and sharp. The destabilizing small-delay sequence makes this a hard negative test: no certificate may hold for mu1 >= mu0 at any small tau.
- **Feasibility:** The 1-D interior version is expressible today (2-D state for the delay). The 2-D boundary version needs a 2-D transport state on an edge, i.e. an L2[s1,s3] state on the edge s2 = 0, and is expressible in principle. The 2-D interior version needs a 3-D state.
- **Citations:**
  - S. Nicaise, C. Pignotti, "Stability and instability results of the wave equation with a delay term in the boundary or internal feedbacks," SIAM J. Control Optim. 45(5):1561-1585, 2006. DOI 10.1137/060648891 [unverified DOI].
  - R. Datko, J. Lagnese, M. P. Polis, "An example on the effect of time delays in boundary feedback stabilization of wave equations," SIAM J. Control Optim. 24:152-156, 1986. DOI 10.1137/0324007 [unverified DOI].

## 8. Delayed heat equation, 1-D (the S18 threshold) and its 2-D extension (Fridman, Orlov 2009)

- **1-D equations [read, Example 1]:**
  - z_t = z_xx + r z - z(x, t - tau(t)) on (0,pi), with z(0) = z(pi) = 0.
  - The paper uses the polytope |r| <= 1.9.
  - Bounds obtained: h = 1.025 for d = 0 (constant delay), and h = 1.021 for an unknown tau_dot.
- **Tightness [read, Remark 4]:** The LMIs coincide with those for the first mode, y_dot = (-1 + r) y - y(t - tau). Stability of the first mode is necessary for stability of the PDE, so the reduction is tight.
- **Exact threshold [derived]:** At r = 1.9 the first mode is y_dot = 0.9 y - y(t - tau). For y_dot = a y - b y(t - tau) with |a| < b, the delay margin is tau* = arccos(a/b) / sqrt(b^2 - a^2) (Hayes 1950 [unverified citation]). Here tau* = arccos(0.9) / sqrt(0.19) = 1.0347. This matches the S18 threshold 1.034726, so S18 is this example with r = 1.9.
- **2-D extension [derived]:**
  - z_t = Delta z + r z - z(t - tau) on (0,pi)^2, with Dirichlet conditions.
  - The first eigenvalue of Delta is -2, so the first mode is y_dot = (r - 2) y - y(t - tau).
  - At r = 2.9 this is the same scalar equation, so tau* = 1.0347.
  - Higher modes have a more negative coefficient and therefore a larger delay margin, since arccos(a/b)/sqrt(b^2-a^2) increases as a decreases.
- **Spaces:**
  - 1-D version: L2[s1] plus a 2-D transport state (= S18).
  - 2-D version: an L2[s1,s2] state plus a 3-D transport state on (0,pi)^2 x [0,tau], whose inflow is the 2-D state itself.
- **Gaps and cases addressed:** (b2) a 2-D delayed PDE; (b4) and (a) a 3-D state with a known threshold that is exactly the S18 number, so the 2-D/3-D result can be checked against the existing 1-D/2-D one.
- **Feasibility:** The 3-D version needs a 3-D state coupled to a 2-D state in-domain, where the 3-D state at theta = 0 equals the 2-D state. Expressible if the 3-D converter accepts a 2-D state as the datum on a face; otherwise missing.
- **Citation:** E. Fridman, Y. Orlov, "Exponential stability of linear distributed parameter systems with time-varying delays," Automatica 45:194-201, 2009. DOI 10.1016/j.automatica.2008.06.006.

## 9. Delayed damped wave, 1-D (Fridman, Orlov 2009, Example 2)

- **Equations [read]:**
  - z_tt = 0.1 z_xx - 2 z_t + u on (0,pi), with Dirichlet ends.
  - The paper adds the delayed feedback u = -z(x, t - tau(t)).
  - Result: exponential stability with decay rate delta = 0.8 for 0 <= tau(t) <= 0.31. This is a sufficient bound obtained at a fixed decay rate.
- **Tightness [read, Remark 7]:** The LMIs of Theorems 5 and 6 coincide with those for the first-mode delay ODE, y_ddot + mu0 y_dot + mu1 y_dot(t - tau) + (a k^2 + a0) y + a1 y(t - tau) = 0 with k = 1. The paper calls the results tight.
- **Exact threshold [derived]:** Here the first mode is y_ddot + 2 y_dot + 0.1 y + y(t - tau) = 0. Its delay margin (without the decay-rate constraint) is computable from the characteristic quasi-polynomial. The paper does not report this number.
- **Spaces:** L2[s1] (first-order wave) plus a 2-D transport state.
- **Case addressed:** (b2), the wave counterpart of S18.
- **Feasibility:** expressible today.
- **Citation:** as item 8.

## 10. Heat equation with a delayed Robin boundary term (Nicaise, Valein, Fridman 2009)

- **Equations [read, system (1)]:**
  - u_t = a u_xx on (0,pi), with u(0,t) = 0.
  - u_x(pi,t) = -mu0 u(pi,t) - mu1 u(pi, t - tau(t)).
  - The paper rewrites the delay as z(rho,t) = u(pi, t - tau(t) rho) on rho in (0,1), with tau z_t + (1 - tau_dot rho) z_rho = 0 and z(0,t) = u(pi,t) [read, system (6)].
- **Threshold:**
  - The paper gives explicit sufficient conditions for exponential stability and the decay rate [abstract].
  - The exact constant-delay margin is not given. It follows from the characteristic equation sqrt(lambda/a) cosh(pi sqrt(lambda/a)) + (mu0 + mu1 e^{-lambda tau}) sinh(pi sqrt(lambda/a)) = 0 [derived, unverified].
- **Spaces:** L2[s1] (heat) plus L2[s2] (transport), coupled at a point.
- **Coupling mechanism:** a point trace of the heat state is the transport inflow, and the outflow enters the Robin condition. This is a boundary-point delay, which the current suite does not have.
- **Case addressed:** (b2).
- **Feasibility:** expressible today for a constant delay. A time-varying delay is not expressible.
- **Citation:** S. Nicaise, J. Valein, E. Fridman, "Stability of the heat and of the wave equations with boundary time-varying delays," Discrete Contin. Dyn. Syst. Ser. S 2(3):559-581, 2009. DOI 10.3934/dcdss.2009.2.559.

## 11. Exact H-infinity gain for self-adjoint strictly negative generators (Lidstrom, Rantzer, Morris 2016)

- **Result [read, Theorem 1]:**
  - The system is z_dot = A z + B u + H d.
  - A is closed, densely defined, self-adjoint and strictly negative. B and H are bounded.
  - The cost output is zeta = (z, u).
  - The minimal closed-loop gain over static state feedback is gamma* = || H* (A^2 + B B*)^{-1} H ||^{1/2}, attained by K_opt = B* A^{-1}.
  - The paper illustrates this on the 1-D heat equation with Dirichlet ends and (H d)(x) = d(t) [read].
- **Application to the gap [derived]:**
  - With B = 0, the open-loop gain from d to z is || A^{-1} H ||.
  - Take A = Delta on (0,1)^2 with Dirichlet conditions. Take a 1-D disturbance d(t,s1) in L2[s1] entering as (H d)(s1,s2) = d(s1). Take the 2-D output z.
  - The sine modes in s1 decouple the problem, so gamma = sup_n || psi_n ||_{L2(0,1)}, where psi_n'' - (n pi)^2 psi_n = 1 and psi_n(0) = psi_n(1) = 0.
  - The supremum is attained at n = 1, with psi_1(s2) = -(1/pi^2) [1 - cosh(pi (s2 - 1/2)) / cosh(pi/2)], so gamma has a closed form.
  - With B nonzero (for example a distributed 2-D input), gamma* gives an exact closed-loop H-infinity value for a synthesis test.
- **Spaces:** a 2-D state, a 1-D disturbance (L2[s1]), and a 2-D output.
- **Coupling mechanism:** a 1-D input entering a 2-D state in-domain, as a multiplier constant in s2.
- **Gap addressed:** (a) an exact gain with a 1-D disturbance and a 2-D output. Analysis and synthesis can both be checked.
- **Feasibility:** expressible today. Boundary inputs are excluded, because the theorem needs a bounded H.
- **Citation:** C. Lidstrom, A. Rantzer, K. Morris, "H-infinity optimal control for infinite-dimensional systems with strictly negative generator," Proc. 55th IEEE CDC, pp. 5275-5280, 2016. DOI 10.1109/CDC.2016.7799077.

## 12. Heat equation with a dynamic boundary condition of reactive-diffusive (Laplace-Beltrami) type

- **Equations [abstract and unverified standard form]:**
  - u_t = Delta u in Omega.
  - On the boundary Gamma: u_t = sigma Delta_Gamma u - d u / d nu (+ beta u).
  - Goldstein 2006 derives general Wentzell conditions for the heat and wave equations from physical principles [abstract]. Whether that paper uses the thin-highly-conductive-layer limit was not confirmed.
- **Box version (proposal) [derived]:**
  - Unit square. The dynamic condition holds on the edge s2 = 1 only, and is carried by a 1-D state v(s1) = u(s1,1) with v_t = sigma v_{s1 s1} - u_{s2}(s1,1) + beta v.
  - At the corners: v_{s1} = 0 at s1 = 0,1, with Neumann conditions in s1 for u.
  - The other edges are Dirichlet.
  - The cosine modes in s1 decouple the problem. Each mode reduces to a 1-D heat equation with a Wentzell (dynamic) condition whose eigenvalues solve a transcendental equation, so a threshold in beta or in a bulk reaction is semi-analytic.
- **Spaces:** L2[s1,s2] plus L2[s1] on an edge.
- **Coupling mechanism:**
  - Two-way.
  - The 1-D state is the Dirichlet datum of the 2-D state, as in S3.
  - The normal derivative of the 2-D state drives the 1-D state, which also carries its own s1-diffusion. This is S4 with tangential diffusion.
- **Case addressed:** (b3) and the bulk-surface gap (a). This is the closest box analogue of the thin-layer and Laplace-Beltrami models.
- **Threshold:** None is published for a box. No paper with explicit eigenvalues on a rectangle was found.
- **Feasibility:** expressible today, since the edge state is second order in s1 with Neumann ends.
- **Citations:**
  - G. R. Goldstein, "Derivation and physical interpretation of general boundary conditions," Adv. Differential Equations 11(4):457-480, 2006 (Project Euclid ade/1355867704; one repository lists pages 419-456).
  - J. L. Vazquez, E. Vitillaro, "Heat equation with dynamical boundary conditions of reactive-diffusive type," J. Differential Equations 250:2143-2161, 2011. DOI 10.1016/j.jde.2010.12.012 [unverified].

## 13. Road-field model: a 1-D line of fast diffusion on the edge of a 2-D field (Berestycki, Roquejoffre, Rossi 2013)

- **Equations [read, eq. (1.1)]:**
  - Road: u_t - D u_xx = nu v(x,0,t) - mu u.
  - Field: v_t - d Delta v = f(v) in y > 0.
  - Exchange: -d v_y(x,0,t) = mu u - nu v(x,0,t).
  - The linearization at 0 replaces f(v) by f'(0) v.
- **Known threshold [read]:** The asymptotic spreading speed on the unbounded domain equals the KPP speed 2 sqrt(d f'(0)) iff D <= 2d, and exceeds it iff D > 2d. This is a speed threshold, not a stability threshold.
- **Box version (proposal) [derived]:**
  - Truncate to (0,L) x (0,1). Put the road on the edge s2 = 0.
  - Use Dirichlet or Neumann conditions elsewhere, and Dirichlet or Neumann conditions at the road ends.
  - Persistence versus extinction of the linearization is decided by the sign of a principal eigenvalue. Berestycki, Ducasse and Rossi define it as the limit of eigenvalues on truncated domains (arXiv:1810.13180) [abstract].
  - On a box with constant coefficients, sine or cosine modes in s1 decouple the problem, so the threshold value of f'(0) solves a 2 x 2 transcendental system.
- **Spaces:** L2[s1,s2] plus L2[s1] on an edge.
- **Coupling mechanism:**
  - Two-way.
  - The 1-D trace v(x,0) enters the road equation in-domain.
  - The road state enters a Robin condition of the 2-D state on the edge.
  - This differs from S3 (Dirichlet datum) and S4 (Wentzell): the exchange is a Robin datum.
- **Cases addressed:** (b3) edge cascades; bulk-surface gap (a).
- **Feasibility:** expressible today.
- **Citations:**
  - H. Berestycki, J.-M. Roquejoffre, L. Rossi, "The influence of a line with fast diffusion on Fisher-KPP propagation," J. Math. Biol. 66:743-766, 2013. arXiv:1210.3721. The journal DOI 10.1007/s00285-012-0508-1 is [unverified].
  - H. Berestycki, R. Ducasse, L. Rossi, "Generalized principal eigenvalues for heterogeneous road-field systems," arXiv:1810.13180.

## 14. Transport-diffusion (hyperbolic-parabolic) couplings and non-coercive Lyapunov functionals (relatives of S17)

- **Ghousein and Witrant 2020 [abstract]:**
  - A diffusion PDE is coupled to a transport PDE in-domain, through an advection source term driven by the transport state and a Volterra integral source term driven by the parabolic state, and at the boundary.
  - Backstepping gives L2-exponential stability of the closed loop. Exact equations were not retrieved (HAL access blocked).
- **Lhachemi, Prieur and Trelat 2026 [abstract]:**
  - A 1-D reaction-diffusion and wave cascade with internal coupling.
  - The open-loop generator is Riesz-spectral under explicit conditions.
  - The paper notes that backstepping for hyperbolic-parabolic systems has mainly treated the case where the hyperbolic part is a transport equation.
- **Mironchenko and Wirth 2019 [abstract]:** converse theorems for non-coercive Lyapunov functions on Banach spaces, with examples that show the assumptions are necessary. This is the theoretical counterpart of the S17 observation that only the non-coercive PIE-to-PDE form certifies.
- **Jacob, Mironchenko, Partington and Wirth [abstract]:** extend non-coercive Lyapunov functions to input-to-state stability (arXiv:1911.01327).
- **No paper found** proves the absence of a coercive quadratic Lyapunov functional for a transport-to-diffusion cascade. This remains open in the search.
- **Spaces:** L2[s1] plus L2[s1] (Ghousein-Witrant, Lhachemi et al.).
- **Mechanisms:** in-domain and boundary coupling between a transport state and a heat state.
- **Case addressed:** (b1).
- **Threshold:** None is known in open loop. These are control designs.
- **Feasibility:** expressible in 1-D. A 2-D analogue of S17 would need the equations from the full text.
- **Citations:**
  - M. Ghousein, E. Witrant, "Backstepping control for a class of coupled hyperbolic-parabolic PDE systems," Proc. ACC 2020, pp. 1600-1605. DOI 10.23919/ACC45564.2020.9147593.
  - H. Lhachemi, C. Prieur, E. Trelat, "Stabilization of a wave-heat cascade system," arXiv:2601.16610, 2026.
  - A. Mironchenko, F. Wirth, "Non-coercive Lyapunov functions for infinite-dimensional systems," J. Differential Equations 266(11):7038-7072, 2019. DOI 10.1016/j.jde.2018.11.026; arXiv:1612.06575.
  - B. Jacob, A. Mironchenko, J. R. Partington, F. Wirth, "Noncoercive Lyapunov functions for input-to-state stability of infinite-dimensional systems," SIAM J. Control Optim. 58(5), 2020, arXiv:1911.01327. The DOI 10.1137/19M1297506 is [unverified].

## 15. Elastic-thermoelastic transmission problem, 1-D (Marzocchi, Munoz Rivera, Naso 2002)

- **Equations [unverified standard form]:**
  - Elastic part on (0,L0): u_tt = a u_xx.
  - Thermoelastic part on (L0,L): v_tt = b v_xx - alpha theta_x and theta_t = k theta_xx - beta v_xt.
  - Transmission at L0: continuity of displacement and stress.
- **Known result [secondary and abstract]:**
  - The energy decays exponentially for this transmission problem.
  - The 2003 IMA paper shows exponential decay of symmetric solutions however small the thermoelastic part is.
  - The second-sound (Cattaneo) variant loses exponential stability (Fernandez Sare, Munoz Rivera, Racke) [abstract].
- **Spaces:** L2[s1] states only, after the reflection.
- **Coupling mechanism:** a point interface between a hyperbolic state and a hyperbolic-parabolic system.
- **Gaps addressed:** (a) thermoelastic coupling and parabolic-hyperbolic coupling. The Cattaneo variant gives a matched positive/negative pair.
- **Threshold:** None (qualitative).
- **Feasibility:** expressible in 1-D. A 2-D strip version on two rectangles has not been found.
- **Citation:** A. Marzocchi, J. E. Munoz Rivera, M. G. Naso, "Asymptotic behaviour and exponential stability for a transmission problem in thermoelasticity," Math. Methods Appl. Sci. 25(11):955-980, 2002. DOI 10.1002/mma.323.

## 16. Age-structured population with spatial diffusion; age-structured chemostat

- **Model [unverified standard linear form]:**
  - p_t + p_a = d p_xx - mu(a) p on (a,x) in (0,A) x (0,L), with Dirichlet or Neumann conditions in x.
  - Renewal condition: p(t,0,x) = int_0^A beta(a) p(t,a,x) da.
- **Threshold [derived for constant d and Dirichlet conditions]:**
  - The sine modes in x decouple the problem.
  - The first mode decays iff R_eff = int_0^A beta(a) exp(-int_0^a mu) exp(-d (pi/L)^2 a) da < 1.
  - Without diffusion this reduces to the net reproduction number R0 = int beta exp(-int mu) < 1, which is the standard McKendrick threshold.
- **Kang 2022 [abstract]:** develops the principal spectral theory for age-structured models with random diffusion and shows that a zero principal eigenvalue is the critical threshold.
- **Karafyllis and Krstic 2017 [abstract]:** treat an age-structured chemostat (a transport state plus a 0-D state), but the model is nonlinear. A linearization would give a 0-D state plus a transport state with an integral renewal condition.
- **Spaces:** a 2-D state on (a,x): transport in s1 = a, diffusion in s2 = x. Optionally an R^q state (a resource) entering the renewal or mortality.
- **Coupling mechanism:** an integral boundary condition (renewal) on the edge a = 0 of a mixed transport-diffusion 2-D state. This is a relative of S17 (first order in one direction, second order in the other).
- **Gap addressed:** (a) age-structured models.
- **Feasibility:** expressible today with polynomial beta and mu (integral boundary condition, mixed-order 2-D state). Exponential survival functions need polynomial surrogates.
- **Citations:**
  - H. Kang, "Effects of diffusion on the principal eigenvalue for age-structured models with random diffusion," Proc. Roy. Soc. Edinburgh A 152(1):258-280, 2022. DOI 10.1017/prm.2021.10.
  - I. Karafyllis, M. Krstic, "Stability of integral delay equations and stabilization of age-structured models," ESAIM: COCV 23(4):1667-1714, 2017. DOI 10.1051/cocv/2016069; arXiv:1609.09437.

## 17. Heat-wave interaction in two and three dimensions (Avalos, Lasiecka, Triggiani 2016)

- **Model [abstract]:**
  - The Stokes part of a fluid-structure model is replaced by an n-dimensional heat equation, giving a vector heat equation coupled to a vector wave equation across an interface.
  - The optimal rational decay rate is proved.
  - The exact rate and the geometry (nested domains, not boxes) were not read [unverified].
- **Spaces:** 2-D or 3-D states.
- **Mechanism:** an interface between two states.
- **Gap addressed:** (a), the 3-D parabolic-hyperbolic case.
- **Feasibility:** The geometry has a structure immersed in a fluid, which is not a box. The rectangular case is covered by item 1 instead.
- **Citation:** G. Avalos, I. Lasiecka, R. Triggiani, "Heat-wave interaction in 2-3 dimensions: optimal rational decay rate," J. Math. Anal. Appl. 437(2):782-815, 2016. DOI 10.1016/j.jmaa.2015.12.045 [unverified DOI].

---

## Not found within the budget

- A published 3-D coupled example on a brick with an explicit stability margin, other than the derived item 8.
- A bulk-surface reaction-diffusion model on a rectangle. Madzvamuse, Chung and Venkataraman 2015 (Proc. R. Soc. A 471:20140546) and Paquin-Lefebvre, Nagata and Ward (SIAM J. Appl. Dyn. Syst., DOI 10.1137/18M1213737) use a disk or sphere.
- A plate-beam (2-D plate with a 1-D beam on an edge) model with a published threshold.
- A paper proving that no coercive quadratic Lyapunov functional exists for a transport-to-diffusion cascade.

## Summary table

| # | Candidate | Spaces | Mechanism | Threshold known | Gap / case | Feasibility |
|---|---|---|---|---|---|---|
| 1 | Wave-heat on two rectangles (BPS 2019) | 2-D x2 (3 after rewrite) | interface traces + normal derivatives on shared edge | qualitative: o(t^{-2/3}), not exp. stable [read] | (a) 2-D x 2-D, 2-D wave, parab.-hyperb. | needs BC coupling two 2-D states on one edge (reflection) |
| 2 | Wave-heat 1-D (BPS 2016; ZZ 2004) | L2[s1] x3 | point interface | qualitative: o(t^{-4}); exact eigen-equation [read] | (a) parab.-hyperb. | today |
| 3 | Wave with acoustic BC (Beale 1976) | 2-D + L2[s1] edge | two-way edge: normal-deriv datum + trace drive | qualitative: not uniformly stable [abstract] | (a) 2-D wave + 1-D; S4/S25 relative | today (first-order wave) |
| 4 | Linear Turing system | 2-D x2 | in-domain multiplier | yes, analytic [derived from standard theory] | (a) 2-D x 2-D in-domain | today |
| 5 | Two-layer heat, rectangles | 2-D x2 | interface flux + Robin jump | semi-analytic [derived] | (a) 2-D x 2-D, piecewise coeffs | needs 2-D/2-D edge BC |
| 6 | Wave, partially delayed boundary control (XYL 2006) | L2[s1] + L2[s2] | point trace -> transport -> Neumann | yes, mu* = 1/2, analytic [abstract] | (b2) | today |
| 7 | Wave with delayed damping (NP 2006; DLP 1986) | 1-D + 2-D transport; 2-D + 3-D | interior or edge delay as transport | yes, mu1 < mu0 sharp [read via Fridman-Orlov] | (b2), (b4) | 1-D today; 2-D needs 3-D/edge transport |
| 8 | Delayed heat, 1-D (= S18) and 2-D extension (FO 2009) | 2-D + 3-D transport | in-domain delay as transport on added dim | yes, tau* = 1.0347 [read 1-D tightness; 2-D derived] | (b2), (b4), (a) 3-D | needs 3-D state with 2-D face datum |
| 9 | Delayed damped wave (FO 2009 Ex. 2) | L2[s1] + 2-D transport | in-domain delay | sufficient 0.31 at decay 0.8 [read]; exact derivable | (b2) | today |
| 10 | Heat with delayed Robin BC (NVF 2009) | L2[s1] + L2[s2] | point trace -> transport -> Robin | sufficient only [abstract]; exact derivable | (b2) boundary-point delay | today (constant delay) |
| 11 | Exact H-inf, self-adjoint negative A (LRM 2016) | 2-D + 1-D input + 2-D output | 1-D disturbance in-domain | yes, closed form [read thm; 2-D derived] | (a) exact gain, 1-D dist., 2-D output | today |
| 12 | Heat with Laplace-Beltrami dynamic BC on one edge | 2-D + L2[s1] edge | two-way: Dirichlet datum + flux drive, tangential diffusion | not published for boxes; semi-analytic [derived] | (b3), (a) bulk-surface | today |
| 13 | Road-field (BRR 2013), boxed | 2-D + L2[s1] edge | two-way Robin exchange | speed threshold D = 2d [read]; box eigen-threshold derived | (b3), (a) bulk-surface | today |
| 14 | Transport-diffusion couplings; non-coercive LFs | L2[s1] x2 | in-domain + boundary | none open-loop | (b1) S17 | 1-D today |
| 15 | Elastic-thermoelastic transmission (MMN 2002) | L2[s1] x several | point interface | qualitative: exp. stable; Cattaneo variant not [abstract] | (a) thermoelastic, parab.-hyperb. | today |
| 16 | Age-structured with diffusion; chemostat | 2-D (transport x diffusion) [+ R^q] | integral renewal BC on edge a = 0 | yes for separable case [derived]; Kang 2022 eigen-threshold [abstract] | (a) age-structured; S17 relative | today (polynomial beta, mu) |
| 17 | Heat-wave 2-3 D (ALT 2016) | 2-D/3-D | interface | rational rate [abstract] | (a) 3-D | non-box geometry; not expressible |

Priority for an addition to the suite:
- Items 1 and 2: a parabolic-hyperbolic interface with a rigorous negative result.
- Items 4 and 11: exact values in the 2-D/2-D and gain gaps.
- Item 8, 2-D version: a 3-D state that reproduces the S18 number.
- Items 12 and 13: two-way edge couplings for case (b3).
- Item 6: a sharp delay threshold.
