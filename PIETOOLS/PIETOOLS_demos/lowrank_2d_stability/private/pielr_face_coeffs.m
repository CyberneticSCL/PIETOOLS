function S = pielr_face_coeffs(V,q,P)                                       % CC, 09/23/2026
% PIELR_FACE_COEFFS  Face coefficients S_i = V_i' X_i V_i in ORIGINAL-b units.
%
% q arrives in NORMALISED-b units (bm_setup divides b by nb0), so the Gram
% block is scaled back here.  Everything pielr_solve RETURNS is in original
% units and the conversion is owned in one place -- getting it wrong reports a
% perfectly good certificate at residual 1.
%
% Same arithmetic as pielr_certify's inline version; separated only because
% pielr_solve is not that file.
S = cell(1,numel(V));
if isempty(V) || isempty(q), return, end
for i = 1:numel(V)
    N  = round(sqrt(numel(P.rows{i})));
    Xi = reshape(q(P.rows{i}),N,N)*P.nb0;
    Si = V{i}'*((Xi+Xi')/2)*V{i};
    S{i} = (Si+Si')/2;
end
end
