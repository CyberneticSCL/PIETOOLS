function val = pi_kernel_pairing(Kcells,xfun,yfun,vars,dom)
% VAL = PI_KERNEL_PAIRING(KCELLS,XFUN,YFUN,VARS,DOM) evaluates the pairing
%
%       <y,Px> = int_s y(s)' (Px)(s) ds
%
% in closed form, for the PI operator on L_2[VARS] whose kernels are given
% as polynomials,
%
%       (Px)(s) = sum_gamma int_{s'} K_gamma(s,s') I_gamma(s-s') x(s') ds'.
%
% INPUTS
% - Kcells: cell array of 3^n3 matrix-valued 'polynomial' kernels, indexed as
%           the parameter cell array of an sopvar object, i.e. cell k holds
%           the kernel for the multi-index PI_GAMMA_INDEX(k,n3). Each kernel
%           is expressed in the variables VARS and their dummy counterparts,
%           named by appending '_dum';
% - xfun:   m x 1 'polynomial' test function in VARS;
% - yfun:   m x 1 'polynomial' test function in VARS;
% - vars:   1 x n3 'cellstr' object naming the spatial variables;
% - dom:    n3 x 2 array specifying the domain of each variable;
%
% OUTPUTS
% - val:    the value of <y,Px>, of type 'double';
%
% Since x, y and every kernel are polynomial, the integrals are evaluated
% exactly by 'int', so the result is not subject to quadrature error.
%
% MP, 08/22/2026: Initial coding

n3 = numel(vars);
vars_dum = strcat(vars,'_dum');

if numel(Kcells)~=3^n3
    error("Expected 3^n3 kernels for n3 spatial variables.")
end

% Express the test function in the dummy variables.
xdum = xfun;
for k = 1:n3
    xdum = subs(xdum,polynomial(vars(k)),polynomial(vars_dum(k)));
end

Px = polynomial(zeros(size(yfun,1),1));
for k = 1:numel(Kcells)
    gam = pi_gamma_index(k,n3);
    Px = Px + apply_indicator(Kcells{k}*xdum,vars_dum,vars,dom,gam);
end

expr = yfun.'*Px;
for k = 1:n3
    expr = int(expr,polynomial(vars(k)),dom(k,1),dom(k,2));
end
val = double(expr);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = apply_indicator(expr,vars_in,vars_out,dom,gam)
% Integrate expr(s,s') against I_gamma(s-s') over s', in closed form.

out = expr;
for k = 1:numel(vars_in)
    spk = polynomial(vars_in(k));
    sk = polynomial(vars_out(k));
    switch gam(k)
        case 1      % delta(s_k-s'_k)
            out = subs(out,spk,sk);
        case 2      % I(s_k-s'_k), i.e. s'_k <= s_k
            out = int(out,spk,dom(k,1),sk);
        case 3      % I(s'_k-s_k), i.e. s'_k >= s_k
            out = int(out,spk,sk,dom(k,2));
    end
end

end
