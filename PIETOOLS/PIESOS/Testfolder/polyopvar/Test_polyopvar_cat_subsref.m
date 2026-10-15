% Test_polyopvar_cat_subsref.m -- verify vertcat and (row-only) '()'
% subsref for 'polyopvar' objects via apply_poly. Oracle: the same
% operations on the numeric polynomials, compared with a relative
% tolerance.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - '()' with an explicit column index (v(2,1)) is not supported:
%   @polyopvar/subsref.m forwards [row,col] to the 'tensopmat' subsref,
%   whose columns are one block per monomial, not per matrix column.
%   Row-only v(2,:) is not affected. The test checks that the column
%   index errors with the message 'number of monomials must match', so
%   it fails once column indices are supported.
% - Tolerance is relative, tol*max(1,|reference|).

tol = 1e-8;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);

xv = {rand_poly([1,1],{'s1','s2'},2,5)};
yv = {rand_poly([1,1],{'s1','s2'},2,5)};

% vertcat
v = [x;y];
if ~isequal(size(v),[2,1])
    error("Test_polyopvar_cat_subsref failed: vertcat size mismatch.")
end
Vv = apply_poly(v,{xv{1},yv{1}});
refv = [xv{1};yv{1}];
errv = cleanpoly(Vv-refv,1e-10);
if max(max(abs(errv.C)))>tol*max(1,max(max(abs(refv.C))))
    error("Test_polyopvar_cat_subsref failed: vertcat value mismatch.")
end

% '()' subsref, row-only
v2 = v(2,:);
Vv2 = apply_poly(v2,{xv{1},yv{1}});
errv2 = cleanpoly(Vv2-yv{1},1e-10);
if max(max(abs(errv2.C)))>tol*max(1,max(max(abs(yv{1}.C))))
    error("Test_polyopvar_cat_subsref failed: row-subsref value mismatch.")
end

% '()' subsref, genuine COLUMN restriction: must error (see NOTES)
f = x + x.*x;
threw = false;
try
    f1 = f(1,1); %#ok<NASGU>
catch ME
    threw = true;
    if ~contains(ME.message,'number of monomials must match')
        error("Test_polyopvar_cat_subsref failed: column-restricted '()' errored with unexpected message: %s",ME.message)
    end
end
if ~threw
    error("Test_polyopvar_cat_subsref failed: column-restricted '()' should have errored (still-open gap, see header note).")
end

disp('Test_polyopvar_cat_subsref passed.')
