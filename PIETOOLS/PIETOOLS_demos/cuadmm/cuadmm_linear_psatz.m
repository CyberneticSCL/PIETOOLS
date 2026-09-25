function S = cuadmm_linear_psatz(S)                                        % CC, 09/25/2026
% S = cuadmm_linear_psatz(S) -- enable the four linear box generators
% eq_use_psatz = [3;4;5;6] (poslpivar_2d, 0cf1add1) on a 2-D settings struct
% returned by cuadmm_settings(tier,2).
%
% A separate file because MATLAB exposes only the first function in a file; a
% local helper inside cuadmm_settings.m would be unreachable by callers.
%
% OFF by default in every cuADMM tier. With a searched P the generators take
% Mosek from failing to certifying at rel_b 1e-08 up to 0.90 lam*, but they add
% +30% to m and 6.6x to nnz(At) (557k -> 3.66M), and nnz(At) is what cuADMM pays
% for on every iteration. cuADMM certified the psatz-OFF 2-D program at
% 0.50 lam* (verified rel_b 2.0e-04), so they are not needed for a certificate
% there. Whether the extra reach is worth the per-iteration cost under cuADMM
% has NOT been measured.
%
% Two traps, both measured:
%  - shipped settings carry only TWO eq_opts_psatz entries (eq_use_psatz is
%    [0;0]), so entries 3 and 4 must be CLONED from an existing one; assigning
%    just .psatz makes a struct with no `exclude` field and poslpivar_2d dies.
%  - each psatz block needs the FULL eq_deg; eq_deg-1 destroys the certificate
%    (reach 0.9999 -> 0.0000, measured).
if ~isfield(S,'settings_2d')
    error('cuadmm_linear_psatz:dim','S has no settings_2d; call cuadmm_settings(tier,2)');
end
pz = [3;4;5;6];
S.settings_2d.eq_use_psatz = pz;
otmpl = S.settings_2d.eq_opts_psatz{1};
for j = 1:numel(pz)
    o = otmpl;  o.psatz = pz(j);
    S.settings_2d.eq_opts_psatz{j} = o;
    S.settings_2d.eq_deg_psatz{j}  = S.settings_2d.eq_deg;
end
S.cuadmm.linear_psatz = true;
end
