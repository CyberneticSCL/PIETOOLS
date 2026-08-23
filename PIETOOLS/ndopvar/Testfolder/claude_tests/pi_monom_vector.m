function Z = pi_monom_vector(Zc,varnames)
% Z = PI_MONOM_VECTOR(ZC,VARNAMES) builds the monomial vector
%
%       Z(x) = kron(Zc{1}(x_1),...,Zc{N}(x_N))
%
% as a 'polynomial' column vector, where Zc{i} is the vector of degrees of
% variable VARNAMES{i}. This is the monomial convention used by the sopvar
% and sdopvar classes for the bases ZL and ZR.
%
% MP, 08/22/2026: Initial coding

degmat = zeros(1,0);
for i = 1:numel(Zc)
    Zi = reshape(Zc{i},[],1);
    degmat = [kron(degmat,ones(numel(Zi),1)), kron(ones(size(degmat,1),1),Zi)];  %#ok<AGROW>
end

if isempty(varnames)
    Z = polynomial(1);
    return
end

T = size(degmat,1);
Z = polynomial(speye(T),degmat,varnames(:),[T,1]);

end
