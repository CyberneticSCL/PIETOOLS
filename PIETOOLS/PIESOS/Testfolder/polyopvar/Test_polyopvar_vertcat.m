% Test_polyopvar_vertcat.m -- verify vertcat for 'polyopvar' objects via apply_poly.
% Oracle: the same operation on the numeric polynomials, compared with a relative
% tolerance tol*max(1,|reference|).

% DJ (CC), 10/15/2026: Initial coding;
% DJ (CC), 10/23/2026: Split from Test_polyopvar_cat_subsref ('()' is in
%   Test_polyopvar_subsref);

tol = 1e-8;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);

xv = {rand_poly([1,1],{'s1','s2'},2,5)};
yv = {rand_poly([1,1],{'s1','s2'},2,5)};

v = [x;y];
if ~isequal(size(v),[2,1])
    error("Test_polyopvar_vertcat failed: size mismatch.")
end
Vv = apply_poly(v,{xv{1},yv{1}});
refv = [xv{1};yv{1}];
errv = cleanpoly(Vv-refv,1e-10);
if max(max(abs(errv.C)))>tol*max(1,max(max(abs(refv.C))))
    error("Test_polyopvar_vertcat failed: value mismatch.")
end

disp('Test_polyopvar_vertcat passed.')