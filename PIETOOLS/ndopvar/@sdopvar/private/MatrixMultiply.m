function params_new = MatrixMultiply(params,L,R)
% params = MatrixMultiply(params,L,R) multiply a sdopvar object by constant matrices
% Cnew=L*C*R, using the identity vec(Cnew) = kron(R',L) * vec(C)
% INPUTS
% -params: parameters (A, Bt) of a sdopvar object C(d), vec(C(d))=A+Bt*d.
% -L:  constant matrix left multiplying P
% -R: constant matrix right multiplying P
%
% OUTPUTS
% - params: parameters of a sdopvar object Cnew(d)=L*C(d)*R
%   vec(Cnew(d))=Anew +Btnew*d, the stored parameters are Anew, Bnew.
K = kron(R', L);
% apply multiplication to each kernel component of sdopvar
for ii = 1:numel(params.A)
    params_new.A{ii} = K * params.A{ii};
    params_new.Bt{ii} = K * params.Bt{ii};
end

end
