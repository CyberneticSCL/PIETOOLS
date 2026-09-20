% DEMO_POINCARE_2D  Certify the 2-D Poincare inequality at low rank.
%
% The inequality:  for u on [0,1]^2 with u(0,.) = u(.,0) = 0,
%     int int a(s1,s2) * u_{s1 s2}^2  >=  c * int int u^2 .
% For a = 1 the sharp constant is c* = (pi^2/4)^2 = 6.088  (extremal
% u = sin(pi s1/2) sin(pi s2/2); confirmed numerically to 9 digits by two
% independent routes during this package's validation).
%
% WHY THIS EXAMPLE, next to the heat demo.  Two reasons:
%   * it is NOT a stability problem -- it exercises pielr_certify_pos, bare
%     positivity of a GIVEN operator, where the certificate is FORCED to
%     realise every cell of the target, including the multiplier cell
%     R22{1,1} = a;
%   * unlike the stability benchmarks, it is genuinely DEGREE-HUNGRY: at the
%     affordable spec below the LPI certifies only a few percent of the true
%     c*, so the degree axis visibly matters here.
%
% DERIVATION of the target operator (validated against semantics on grids,
% with a deliberately-wrong negative control, before any rank was measured).
% Put v = u_{s1 s2}, so u = (Kv)(s1,s2) = int_0^{s1} int_0^{s2} v.  Then
%     int a u_{s1s2}^2 - c int u^2  =  <v, (M_a - c K'K) v>,
%     (K'K v)(s1,s2) = int int v(t1,t2) (1-max(s1,t1)) (1-max(s2,t2)) dt,
% and splitting the kernel by max-branch gives the opvar2d cells below:
% R22{1,1} the multiplier a, R22{2,2}/{2,3}/{3,2}/{3,3} the four integral
% branches (index 2 = lower integral, 3 = upper, per direction).
%
% MEASURED at spec [2 2 1], a = 1 (this package's validation runs): the LPI's
% own certifiable boundary is c_max = 0.2318 (3.8% of c*), and at
% c = 0.5*c_max the minimum rank found is 3 of N = 153 -- 11,781 Gram unknowns
% reduced to 6 -- with the rejected rank 2 failing at ~1e-5, a genuine
% near-obstruction.  Try a non-constant weight afun and watch rank(Q11) track
% the number of squares the weight needs (theory.pdf, Section 7).
%
% CC, 09/19/2026: initial version, porting the validated forcing-test builder;
%                 written as a SCRIPT so cert, Ptgt, ... persist in the
%                 workspace for further tests (maintainer request).

here = fileparts(mfilename('fullpath'));
cd(here);
pielr_path();

% ---- the inequality to certify ----------------------------------------------
afun = @(x,y) 1 + 0*x;      % the weight a(s1,s2); try 1 + x.^2 + y.^2
c    = 0.1159;              % operating point: 0.5 * measured c_max at [2 2 1]

pvar s1 s2 s1_dum s2_dum
Ptgt = opvar2d([],[0 0;0 0;0 0;1 1],[0 1;0 1],[s1;s2],[s1_dum;s2_dum]);
Ptgt.R22{1,1} = polynomial(afun(s1,s2)) + 0*s1;          % the multiplier: FORCED
Ptgt.R22{2,2} = -c*(1-s1)*(1-s2)         + 0*s1;         % t1<s1, t2<s2
Ptgt.R22{2,3} = -c*(1-s1)*(1-s2_dum)     + 0*s1;         % t1<s1, t2>s2
Ptgt.R22{3,2} = -c*(1-s1_dum)*(1-s2)     + 0*s1;         % t1>s1, t2<s2
Ptgt.R22{3,3} = -c*(1-s1_dum)*(1-s2_dum) + 0*s1;         % t1>s1, t2>s2

% ---- certify ------------------------------------------------------------------
cert = pielr_certify_pos(Ptgt);       % defaults: spec [2 2 1], BM discovery
assert(cert.ok,'demo_poincare_2d: did not certify -- see the report above.');

fprintf(['\nCertified  int int a*u_{s1s2}^2 >= %.4f * int int u^2  at Gram ' ...
         'rank %d of %d\n(true c* for a = 1 is (pi^2/4)^2 = %.4f; the gap is ' ...
         'the LPI degree, not the rank).\n'],c,cert.r,cert.Ns(1),(pi^2/4)^2);
fprintf('DEMO_POINCARE_2D_DONE\n');
