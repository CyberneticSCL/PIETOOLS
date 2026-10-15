% Test_polyopvar_times_single_nonzero.m -- @polyopvar/times.m for B with a
% single nonzero in B.Z (e.g. a single state x) and A with several
% monomials, e.g. (x+2*x^2).*x: the value is compared with the explicit
% product of the polynomials in the values of the states.
%
% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-10;
pvar s
polyopvar x(s) y(s)
xv = 2+3*s;     yv = 1-s^2;
ncases = 0;

cases = {   (x+2*(x.*x)).*x,                    (xv+2*xv^2)*xv
            x.*(x+2*(x.*x)),                    xv*(xv+2*xv^2)
            (x+(x.*x)+3*(x.*x.*x)).*x,          (xv+xv^2+3*xv^3)*xv
            (x+y).*x,                           (xv+yv)*xv
            (x+y).*y,                           (xv+yv)*yv
            (x+2*y+x.*y).*y,                    (xv+2*yv+xv*yv)*yv
            (x+y).*(x+y),                       (xv+yv)^2     };
for k=1:size(cases,1)
    Y = apply_poly(cases{k,1},{xv,yv});
    D = Y-cases{k,2};
    if ~isempty(nonzeros(D.coefficient)) && max(abs(nonzeros(D.coefficient)))>tol
        error("Test_polyopvar_times_single_nonzero failed (%d): residual %.3e.",k,max(abs(nonzeros(D.coefficient))))
    elseif isempty(nonzeros(polynomial(cases{k,2}).coefficient))
        error("Test_polyopvar_times_single_nonzero (%d): reference value is zero; test would be vacuous.",k)
    end
    ncases = ncases+1;
end

fprintf('Test_polyopvar_times_single_nonzero passed (%d cases).\n',ncases)
