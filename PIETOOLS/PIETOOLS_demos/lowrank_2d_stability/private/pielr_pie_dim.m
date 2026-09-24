function d = pielr_pie_dim(PIE)                                             % CC, 09/23/2026
% PIELR_PIE_DIM  1 or 2, from a pie_struct's own field or from the class of T.
if isa(PIE,'pie_struct'), d = PIE.dim; return, end
if isstruct(PIE) && isfield(PIE,'dim') && ~isempty(PIE.dim) && isscalar(PIE.dim)
    d = PIE.dim;  return
end
if isa(PIE.T,'opvar2d') || isa(PIE.T,'dopvar2d'), d = 2; else, d = 1; end
end
