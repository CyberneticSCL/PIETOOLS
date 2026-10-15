% Test_polyopvar_subsref_curly.m -- verify @polyopvar/subsref.m 'case {}':
% f{i,j} must return a standalone 'polyopvar' object representing the term
% f.C.ops{i,j} paired with monomial f.Z(j,:), i.e. the (i,j) contribution
% to f(x) = C*Z(x). Summing apply_poly over every populated column j,
% for a fixed row i, must reconstruct apply_poly(f,.) restricted to
% row i exactly.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - f{i,j} must return a 'tensopmat' coefficient (not the bare
%   tensopvar/intop that C.ops{i,j} unwraps to), else out.dims.out /
%   out.depmat.out fail.
% - A (row,col) pair whose own block is empty and which is the sole entry
%   of its extracted row and column has no inferable dimensions:
%   tensopmat/initialize reports "ambiguous". Checked as an expected
%   error, not a defect.

tol = 1e-8;
ntrials = 5;

pvar s1 s2
x = polyopvar('x',[s1,s2],[0,1]);
y = polyopvar('y',[s1,s2],[0,1]);
f = x + x.*x;       % 2 distinct monomial terms: bare x, and x^2
v = [f;y];          % row 1: f's 2 terms; row 2: y's single (different) term

[nr,nc] = size(v.C);
if nc<2
    error("Test_polyopvar_subsref_curly failed: test setup expected >=2 monomial columns.")
end

for trial=1:ntrials

xv = rand_poly([1,1],{'s1','s2'},2,5);
yv = rand_poly([1,1],{'s1','s2'},2,5);
Yfull = apply_poly(v,{xv,yv});

for i=1:nr
    Ysum = 0;
    for j=1:nc
        op_ij = v.C.ops{i,j};
        if isempty(op_ij)
            continue
        end
        vij = v{i,j};
        if ~isa(vij,'polyopvar')
            error("Trial %d: v{%d,%d} did not return a 'polyopvar' object.",trial,i,j)
        elseif ~isequal(vij.Z,v.Z(j,:))
            error("Trial %d: v{%d,%d}.Z does not match v.Z(%d,:).",trial,i,j,j)
        end
        Ysum = Ysum + apply_poly(vij,{xv,yv});
    end
    erri = cleanpoly(Ysum-Yfull(i),1e-10);
    if max(max(abs(erri.C)))>tol*max(1,max(max(abs(Yfull(i).C))))
        error("Trial %d: sum of extracted terms for row %d does not reconstruct the original row.",trial,i)
    end
end

end

%%%% Linear indexing must match the corresponding (row,col) pair
v11 = v{1,1};   vlin1 = v{1};
if ~isequal(v11.Z,vlin1.Z) || ~isequal(v11.C.ops,vlin1.C.ops)
    error("Test_polyopvar_subsref_curly failed: linear index v{1} does not match v{1,1}.")
end

%%%% Out-of-range indices must error
threw = false;
try
    bad = v{1,nc+1}; %#ok<NASGU>
catch
    threw = true;
end
if ~threw
    error("Test_polyopvar_subsref_curly failed: out-of-range column index should have errored.")
end
threw = false;
try
    bad = v{nr+1,1}; %#ok<NASGU>
catch
    threw = true;
end
if ~threw
    error("Test_polyopvar_subsref_curly failed: out-of-range row index should have errored.")
end

%%%% A (row,col) pair whose sole (1x1) block is empty cannot be extracted
empty_pair_found = false;
for i=1:nr
    for j=1:nc
        if isempty(v.C.ops{i,j})
            empty_pair_found = true;
            threw = false;
            try
                bad = v{i,j}; %#ok<NASGU>
            catch
                threw = true;
            end
            if ~threw
                error("Test_polyopvar_subsref_curly failed: extracting an all-empty (%d,%d) block should have errored.",i,j)
            end
        end
    end
end
if ~empty_pair_found
    error("Test_polyopvar_subsref_curly failed: test setup expected at least one empty (row,col) pair.")
end

disp('Test_polyopvar_subsref_curly passed.')
