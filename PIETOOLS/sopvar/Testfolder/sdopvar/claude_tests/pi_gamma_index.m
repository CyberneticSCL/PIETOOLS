function idx = pi_gamma_index(k,n3)
% IDX = PI_GAMMA_INDEX(K,N3) converts the linear index K of a sopvar/sdopvar
% parameter cell into the corresponding multi-index, with entries in
% {1,2,3} denoting respectively a multiplier (delta), a lower integral and
% an upper integral in each of the N3 spatial variables.
%
% The convention matches the 3 x 3 x ... x 3 layout of the parameter cell
% array, in which the first spatial variable is the fastest index.
%
% MP, 08/22/2026: Initial coding

idx = ones(1,n3);
k = k-1;
for i = 1:n3
    idx(i) = mod(k,3)+1;
    k = floor(k/3);
end
if k~=0
    error("Linear index exceeds 3^n3.")
end

end
