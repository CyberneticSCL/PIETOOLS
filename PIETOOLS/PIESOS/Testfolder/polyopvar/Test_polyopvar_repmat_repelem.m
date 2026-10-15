% Test_polyopvar_repmat_repelem.m -- verify repmat.m/repelem.m (row-wise
% tiling and per-row replication) via apply_poly.
% Oracle: repmat/repelem of the numeric vector of state values.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - Column replication (dimension 2) is not tested: both methods reject it,
%   as matdim(2)>1 is not supported in this class (see
%   Test_polyopvar_horzcat.m).

tol = 1e-8;

for trial = 1:20
    rng(5000+trial);

    pvar s1 s2
    x = polyopvar('x',[s1,s2],[0,1]);
    y = polyopvar('y',[s1,s2],[0,1]);
    v = [x;y];

    xv = rand_poly([1,1],{'s1','s2'},2,5);
    yv = rand_poly([1,1],{'s1','s2'},2,5);
    base = [xv;yv];

    % repmat, [n,1]
    n = randi([1,4]);
    r = repmat(v,[n,1]);
    if ~isequal(size(r),[2*n,1])
        error("Test_polyopvar_repmat_repelem failed: repmat size mismatch (trial=%d).",trial)
    end
    Yr = apply_poly(r,{xv,yv});
    expected_r = repmat(base,n,1);
    err_r = cleanpoly(Yr-expected_r,1e-10);
    if max(max(abs(err_r.C)))>tol
        error("Test_polyopvar_repmat_repelem failed: repmat value mismatch (trial=%d).",trial)
    end

    % repmat, scalar form repmat(v,n) means [n,n] -- column reps~=1 must error
    try
        repmat(v,2);
        error("Test_polyopvar_repmat_repelem failed: repmat(v,2) should have errored (trial=%d).",trial)
    catch ME
        if contains(ME.message,'Test_polyopvar_repmat_repelem')
            rethrow(ME);
        end
    end

    % repelem, vector R
    R = randi([0,3],2,1);
    if sum(R)==0
        R(1) = 1; % avoid the (also-rejected) all-zero case
    end
    re = repelem(v,R);
    if ~isequal(size(re),[sum(R),1])
        error("Test_polyopvar_repmat_repelem failed: repelem size mismatch (trial=%d).",trial)
    end
    Yre = apply_poly(re,{xv,yv});
    expected_re = repelem(base,R,1);
    err_re = cleanpoly(Yre-expected_re,1e-10);
    if max(max(abs(err_re.C)))>tol
        error("Test_polyopvar_repmat_repelem failed: repelem value mismatch (trial=%d).",trial)
    end

    % repelem, scalar R
    Rs = randi([1,3]);
    re2 = repelem(v,Rs);
    Yre2 = apply_poly(re2,{xv,yv});
    expected_re2 = repelem(base,Rs,1);
    err_re2 = cleanpoly(Yre2-expected_re2,1e-10);
    if max(max(abs(err_re2.C)))>tol
        error("Test_polyopvar_repmat_repelem failed: repelem scalar value mismatch (trial=%d).",trial)
    end
end

disp('Test_polyopvar_repmat_repelem passed.')
