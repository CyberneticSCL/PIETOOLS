function G = gate1d(prog,H,P,q,do_gal,NG)
% PROVENANCE.  scratchpad/reach1d/gate1d.m verbatim.  It is the 1-D analogue of
% the shipped private/opcheck_2d: residual AND per-block PSD, never a solver
% status flag.  T1 is the test that exists because EACH HALF ALONE has
% manufactured a wrong answer here (defect B4).
% NOTE the R1_GATE environment override below: T0 asserts it is unset, because
% a stray value silently redefines what every other test means.
%
% gate1d(prog,H,P,q) -- THE acceptance test for this study, 1-D.
%
% Core is the banked scratchpad/opcheck.m verbatim (copied into reach1d/): it
% substitutes the Gram vector q -- NORMALISED-b units, as bm_setup/bm_resid
% produce -- back through the real opvar operators, rebuilds
% Dop = T'PA + A'PT + epneg T'PT from the SUBSTITUTED Pop, and reports
%   relP = max|Dop+Deop| / max|Pop|
% together with each Gram block's min eigenvalue in original units.  relP, not
% rel(Dop): dividing by max|Dop| reads 0/0 wherever a candidate drives Dop to
% zero as an operator (measured on the Example 21 beam), while max|Pop| cannot
% collapse because Pop >= eppos2*I by construction.
%
% ACCEPTANCE (never a solver status flag):
%   (a) relP < 1e-6
%   (b) every Gram block PSD:  mineig_i >= -1e-8 * ||Q_i||_F
% and, when do_gal is set, two Gram-free operator confirmations:
%   (c) min eig of the Gauss-Legendre Galerkin matrix of Pop  > 0
%   (d) min eig of the Galerkin matrix of -Dop  > -1e-6*max|Pop|
% (c)-(d) test the operators themselves rather than the Gram parametrisation,
% so they cannot be satisfied by a Gram artefact; they are the 1-D form of the
% checks ex21_deg.m applies.  They are optional only because they cost a
% second substitution, not because they are weaker.
if nargin<5||isempty(do_gal), do_gal = false; end
if nargin<6||isempty(NG),     NG = 24; end
R = opcheck(prog,H,P,q);
G = R;
G.psd  = all(R.mineig >= -1e-8*max(R.normQ,realmin));
% The acceptance threshold is 1e-6, the value used throughout this campaign.
% R1_GATE overrides it, for the one experiment that asks how much of a
% measured reach is real and how much is scraping the threshold: BM's accepted
% points sit at relP ~ 1e-7 while the IPM's sit at 1e-9..1e-10, so the two
% curves are not equally far inside the gate and a reach read at one threshold
% must be read again at a tighter one before it is believed.
gt = str2double(getenv('R1_GATE'));   if isnan(gt)||gt<=0, gt = 1e-6; end
G.gate_thresh = gt;
G.cert = (R.relP < gt) && G.psd;
G.ok   = G.cert;                  % restrict_solve_1d's caller convention
G.gal  = false;  G.ePmin = NaN;  G.eDmin = NaN;
if do_gal
    pg = prog;  pg.solinfo.info = struct('verified_by','gate1d');
    pg.solinfo.RRx = full(q(:))*P.nb0;
    Pop  = lpigetsol(pg,H.Pop);
    Deop = lpigetsol(pg,H.Deop);                                     %#ok<NASGU>
    Dop  = H.Top'*Pop*H.Aop + H.Aop'*Pop*H.Top + H.st.epneg*H.Top'*Pop*H.Top;
    eP = eig(op2gal(Pop,NG));   eD = eig(op2gal(-Dop,NG));
    G.ePmin = min(eP);   G.eDmin = min(eD);
    G.gal = (G.ePmin > 0) && (G.eDmin > -1e-6*R.maxPop);
    G.cert = G.cert && G.gal;
end
end
