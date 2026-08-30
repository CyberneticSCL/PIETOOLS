function Kcells = pi_sdopvar_kernels(P,dval)
% KCELLS = PI_SDOPVAR_KERNELS(P,DVAL) evaluates the kernels of an 'sdopvar'
% object at the given values of its decision variables, returning them as
% matrix-valued 'polynomial' objects,
%
%       K_gamma(s,s') = (I_m kron ZL(s))' unvec(A_gamma+B_gamma'*d) ...
%                                                       (I_n kron ZR(s')).
%
% INPUTS
% - P:      'sdopvar' object;
% - dval:   numel(P.Zd) x 1 array of values for the decision variables, in
%           the order of P.Zd;
%
% OUTPUTS
% - Kcells: cell array of the same size as P.params.A, holding the kernel of
%           each parameter cell as a P.dims(1) x P.dims(2) 'polynomial' in
%           the output variables P.vars.out and the dummy input variables,
%           the latter named by appending '_dum' to P.vars.in;
%
% MP, 08/22/2026: Initial coding

if ~isa(P,'sdopvar')
    error("Input must be an 'sdopvar' object.")
end
dval = reshape(dval,[],1);
if numel(dval)~=numel(P.Zd)
    error("Number of decision variable values must match numel(P.Zd).")
end

m1 = P.dims(1);     m2 = P.dims(2);
NL = prod(cellfun(@numel,P.ZL));
NR = prod(cellfun(@numel,P.ZR));

ZLp = pi_monom_vector(P.ZL,P.vars.out);
ZRp = pi_monom_vector(P.ZR,strcat(P.vars.in,'_dum'));

Kcells = cell(size(P.params.A));
for k = 1:numel(P.params.A)
    Ck = reshape(full(P.params.A{k} + P.params.B{k}.'*dval),m1*NL,m2*NR);
    Kk = polynomial(zeros(m1,m2));
    for p = 1:m1
        for q = 1:m2
            blk = Ck((p-1)*NL+(1:NL),(q-1)*NR+(1:NR));
            Kk(p,q) = ZLp.'*(blk*ZRp);
        end
    end
    Kcells{k} = Kk;
end

end
