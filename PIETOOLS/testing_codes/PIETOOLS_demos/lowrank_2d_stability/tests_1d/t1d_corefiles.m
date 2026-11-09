function nm = t1d_corefiles()                                               % CC, 09/22/2026
% t1d_corefiles -- the SEVEN dimension-agnostic files this suite exists to test.
%
% These are the Burer-Monteiro core of the shipped low-rank certifier.  Nothing
% in them mentions a spatial dimension: they see only the SeDuMi triple
% (Atf,bf,Ns,Kf) and a rank profile.  Every defect this campaign found lived
% either in one of them or in restrict_solve, which is why a 1-D driver over
% this same core is a test of the shipped code and not of a re-implementation
% -- and why T0's byte-identity assertion is the load-bearing test of the
% suite.  If these diverge from ../private/, every other result here becomes a
% statement about a copy, so T0 fails the whole run rather than reporting one
% failed test among ten.
nm = {'bm_setup','bm_dr','bm_lm2','bm_proj','bm_resid','bm_report','raw_data'};
end
