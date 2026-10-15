% Test_polyopvar_common_basis.m -- verify common_basis.m/common_vars.m
% (used internally by plus/horzcat/vertcat/times) correctly re-express two
% 'polyopvar' objects, built from disjoint sets of spatial variables and
% state names, over a shared monomial basis and shared variable list,
% without changing the value each one represents. Oracle: the value of each
% operand before, by apply_poly on its own states, equals its value after.

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - States have disjoint spatial variables to exercise the variable-union
%   logic; '+' itself cannot be used across differing dependencies (see
%   Test_polyopvar_arith.m).

tol = 1e-8;

for trial = 1:20
    rng(3000+trial);

    pvar s1 s2 s3
    x = polyopvar('x',s1,[0,1]);
    y = polyopvar('y',[s2,s3],[0,1;0,1]);

    A = x.*x;           % depends only on x (hence s1)
    B = y.*y;           % depends only on y (hence s2,s3)

    xv = rand_poly([1,1],{'s1'},2,4);
    yv = rand_poly([1,1],{'s2','s3'},2,4);

    [A2,B2] = common_basis(A,B);

    % shared basis and variable list
    if ~isequal(A2.varname,B2.varname)
        error("Test_polyopvar_common_basis failed: varname not unified (trial=%d).",trial)
    end
    if ~isequal(full(A2.Z),full(B2.Z))
        error("Test_polyopvar_common_basis failed: Z basis not unified (trial=%d).",trial)
    end
    if ~ismember('x',A2.varname) || ~ismember('y',A2.varname)
        error("Test_polyopvar_common_basis failed: unified varname missing a state (trial=%d).",trial)
    end

    % values unchanged (A2 depends only on x, B2 only on y, though over the same 2-state basis)
    [~,ix] = ismember('x',A2.varname);
    [~,iy] = ismember('y',A2.varname);
    xvals = cell(2,1);
    xvals{ix} = xv;
    xvals{iy} = yv;

    YA_orig = apply_poly(A,{xv});
    YA_new  = apply_poly(A2,xvals);
    errA = cleanpoly(YA_orig-YA_new,1e-10);
    if max(max(abs(errA.C)))>tol
        error("Test_polyopvar_common_basis failed: A value changed by common_basis (trial=%d).",trial)
    end

    YB_orig = apply_poly(B,{yv});
    YB_new  = apply_poly(B2,xvals);
    errB = cleanpoly(YB_orig-YB_new,1e-10);
    if max(max(abs(errB.C)))>tol
        error("Test_polyopvar_common_basis failed: B value changed by common_basis (trial=%d).",trial)
    end
end

disp('Test_polyopvar_common_basis passed.')
