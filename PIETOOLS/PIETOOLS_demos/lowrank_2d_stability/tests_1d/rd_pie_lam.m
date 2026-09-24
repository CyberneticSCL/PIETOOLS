function PIE = rd_pie_lam(n,lam)
% PROVENANCE.  scratchpad/reach1d/rd_pie_lam.m verbatim.  The system under
% test: stable iff lam < lam* = pi^2 = 9.8696, so the suite has a KNOWN
% boundary to test either side of -- which is what makes T8's attribution
% check non-vacuous (above lam* no certificate exists, so a miss there must be
% charged to the relaxation, not to the search).
%
% rd_pie_lam(n,lam) -- 1-D reaction-diffusion  phi_t = phi_ss + lam*phi on
% [0,1], Dirichlet, n identical decoupled states.  Identical to the banked
% ladder factory mk_pie('rd1d',n) / rd_pie(n) except that the reaction
% coefficient is an argument instead of the hard-coded 2, so lam can be swept
% up to the analytic limit lam* = pi^2 (eigenvalues lam - j^2 pi^2).
pvar s t
phi = pde_var('state',n,s,[0,1]);
sys = [diff(phi,t,1)==diff(phi,s,2)+lam*phi;
       subs(phi,s,0)==zeros(n,1);
       subs(phi,s,1)==zeros(n,1)];
PIE = convert(sys);
end
