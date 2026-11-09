function f = lb_fp(op)
% Fingerprint of a 1-D opvar: sum, sum of |.| and number of monomials of
% the coefficients of P, Q1, Q2, R0, R1, R2.
parts = {op.P, op.Q1, op.Q2, op.R.R0, op.R.R1, op.R.R2};
f = zeros(1,3*numel(parts));
for k = 1:numel(parts)
    x = parts{k};
    if isa(x,'polynomial')
        c = x.coefficient;   nm = size(x.degmat,1);
    else
        c = x;   nm = 0;
    end
    f(3*k-2:3*k) = [full(sum(c(:))), full(sum(abs(c(:)))), nm];
end
end
