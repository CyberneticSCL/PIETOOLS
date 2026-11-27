function v = d3_field(K,f)
% V = D3_FIELD(K,F), 11/27/2026: field F of the cone struct K, [] if absent.
v = [];
if isstruct(K) && isfield(K,f),     v = K.(f);  end
end
