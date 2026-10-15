% Test_polyopvar_arith.m -- verify plus/minus/uminus/times/mtimes for
% 'polyopvar' objects via apply_poly, using two state variables sharing
% the same spatial dependency. Oracle: the same operation on the numeric
% polynomials, compared with a relative tolerance.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - 'x+y' for states with different spatial dependencies is not tested:
%   it errors in @tensopmat/plus.m (common_basis.m does not form a union
%   depmat for the sum).
% - Tolerance is relative, tol*max(1,|reference|), as output magnitude is
%   not bounded a priori.

tol = 1e-8;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);

xv = {rand_poly([1,1],{'s1','s2'},2,5)};
yv = {rand_poly([1,1],{'s1','s2'},2,5)};

% plus
z = x+y;
Zv = apply_poly(z,{xv{1},yv{1}});
ref = xv{1}+yv{1};
err = cleanpoly(Zv-ref,1e-10);
if max(max(abs(err.C)))>tol*max(1,max(max(abs(ref.C))))
    error("Test_polyopvar_arith failed: plus mismatch.")
end

% minus
z2 = x-y;
Zv2 = apply_poly(z2,{xv{1},yv{1}});
ref2 = xv{1}-yv{1};
err2 = cleanpoly(Zv2-ref2,1e-10);
if max(max(abs(err2.C)))>tol*max(1,max(max(abs(ref2.C))))
    error("Test_polyopvar_arith failed: minus mismatch.")
end

% uminus
z3 = -x;
Zv3 = apply_poly(z3,{xv{1},yv{1}});
err3 = cleanpoly(Zv3+xv{1},1e-10);
if max(max(abs(err3.C)))>tol*max(1,max(max(abs(xv{1}.C))))
    error("Test_polyopvar_arith failed: uminus mismatch.")
end

% times (scalar .* polynomial)
a = 3.5;
z4 = a.*x;
Zv4 = apply_poly(z4,{xv{1},yv{1}});
ref4 = a*xv{1};
err4 = cleanpoly(Zv4-ref4,1e-10);
if max(max(abs(err4.C)))>tol*max(1,max(max(abs(ref4.C))))
    error("Test_polyopvar_arith failed: times mismatch.")
end

% mtimes (scalar * polynomial)
z5 = a*x;
Zv5 = apply_poly(z5,{xv{1},yv{1}});
ref5 = a*xv{1};
err5 = cleanpoly(Zv5-ref5,1e-10);
if max(max(abs(err5.C)))>tol*max(1,max(max(abs(ref5.C))))
    error("Test_polyopvar_arith failed: mtimes mismatch.")
end

disp('Test_polyopvar_arith passed.')
