% Test_polyopvar_dmonomials.m -- verify sosprogramming/dmonomials.m (used
% by several Examples/ files and sosprogramming/piesos_ineq.m to build
% monomial vectors) against direct polynomial evaluation, for varying
% numbers of state variables and degrees. Oracle: the number of monomials
% (nchoosek sum) and the values of the monomials at a random point,
% compared as sorted sets (the order of the monomials in Z is not
% asserted).

% DJ (CC), 10/15/2026: Initial coding;

% NOTES:
% - Only correctness is tested. dmonomials.m builds each monomial by
%   repeated '*' and carries the full monomial basis seen so far, without
%   pruning zero blocks (e.g. 7 basis rows for a 6-entry degree<=2 vector
%   in 2 variables, one with an all-zero coefficient).

tol = 1e-8;

for nvars = 1:3
    for trial = 1:8
        rng(7000*nvars+trial);

        pvar s1 s2
        xs = cell(nvars,1);
        for i=1:nvars
            xs{i} = polyopvar(['x',num2str(i)],[s1,s2],[0,1]);
        end
        vartab = vertcat(xs{:});

        xvals = cell(nvars,1);
        for i=1:nvars
            xvals{i} = rand_poly([1,1],{'s1','s2'},2,4);
        end

        maxdeg = randi([0,3]);
        degs = (0:maxdeg)';
        Z = dmonomials(vartab,degs);

        % expected count = sum_{d=0}^{maxdeg} nchoosek(nvars+d-1,d)
        expected_count = 0;
        for d=0:maxdeg
            expected_count = expected_count + nchoosek(nvars+d-1,d);
        end
        if size(Z,1)~=expected_count
            error("Test_polyopvar_dmonomials failed: monomial count mismatch (nvars=%d,maxdeg=%d,trial=%d).",nvars,maxdeg,trial)
        end

        Y = apply_poly(Z,xvals);

        % values at a random point, compared as sorted sets
        s1val = rand(); s2val = rand();
        xnum = zeros(nvars,1);
        for i=1:nvars
            xnum(i) = double(subs(xvals{i},[s1;s2],[s1val;s2val]));
        end
        Ynum = zeros(expected_count,1);
        for k=1:expected_count
            Ynum(k) = double(subs(Y(k),[s1;s2],[s1val;s2val]));
        end
        % expected values: all compositions of d into nvars nonneg integer parts
        expected_num = [];
        for d=0:maxdeg
            combos = compositions_nonneg(d,nvars);
            for c=1:size(combos,1)
                expected_num = [expected_num; prod(xnum(:)'.^combos(c,:))];
            end
        end
        if numel(expected_num)~=numel(Ynum)
            error("Test_polyopvar_dmonomials failed: expected/actual count mismatch (nvars=%d,maxdeg=%d,trial=%d).",nvars,maxdeg,trial)
        end
        d1 = sort(expected_num);
        d2 = sort(Ynum);
        if max(abs(d1-d2)) > tol*max(1,max(abs(d1)))
            error("Test_polyopvar_dmonomials failed: value mismatch (nvars=%d,maxdeg=%d,trial=%d).",nvars,maxdeg,trial)
        end
    end
end

disp('Test_polyopvar_dmonomials passed.')


function combos = compositions_nonneg(d,n)
% All n-tuples of nonnegative integers summing to d (rows of an array).
if n==1
    combos = d;
    return
end
combos = zeros(0,n);
for k=0:d
    sub = compositions_nonneg(d-k,n-1);
    combos = [combos; [k*ones(size(sub,1),1), sub]];
end
end
