% Test_polyopvar_vector_states.m -- verify 'polyopvar' objects built from
% vector-valued state variables (dims.in(i)>1) are handled correctly by
% plus/times/mtimes/vertcat/combine_terms, numerically via apply_poly --
% not just structurally -- and that objects of 3 or more states with
% mixed (some scalar, some vector) sizes work as well.
%
% DJ (CC), 10/15/2026: Initial coding;
%
% NOTES:
% - The convenience constructor polyopvar(name,vars,dom) sets dims.in(i)=1
%   for every state, so the states here come from rand_polyopvar
%   (Testfolder/polyopvar).
% - initialize infers dims.in(i) from the coefficient operator; scalar is
%   only the fallback when a state never appears as a bare degree-1
%   monomial (Part 2).

tol = 1e-8;

%% Part 1: arithmetic on vector/mixed-size multi-state objects
for trial = 1:20
    rng(11000+trial);

    [x,xvarname,dims_in,vars,dom,depmat_in] = rand_polyopvar();
    p = numel(xvarname);

    xvals = cell(p,1);
    for i=1:p
        xvals{i} = rand_poly([dims_in(i),1],vars(depmat_in(i,:)),2,5);
    end

    % Identity evaluation: x should literally evaluate to the stacked xvals
    Y = apply_poly(x,xvals);
    expected = vertcat(xvals{:});
    err = cleanpoly(Y-expected,1e-10);
    if max(max(abs(err.C)))>tol
        error("Test_polyopvar_vector_states failed: identity eval mismatch (trial=%d,p=%d,dims_in=%s).", ...
            trial,p,mat2str(dims_in'))
    end

    % Elementwise square x.*x of the full (possibly multi-state, possibly
    % mixed-size) stacked vector must stay the same overall size
    x2 = x.*x;
    if ~isequal(size(x2),[sum(dims_in),1])
        error("Test_polyopvar_vector_states failed: size(x.*x) mismatch (trial=%d).",trial)
    end
    Y2 = apply_poly(x2,xvals);
    err2 = cleanpoly(Y2-expected.*expected,1e-10);
    if max(max(abs(err2.C)))>tol
        error("Test_polyopvar_vector_states failed: x.*x value mismatch (trial=%d).",trial)
    end

    % Scalar linear combination (mtimes/plus)
    a = randn; b = randn;
    z3 = a*x + b*x;
    Y3 = apply_poly(z3,xvals);
    err3 = cleanpoly(Y3-(a+b)*expected,1e-10);
    if max(max(abs(err3.C)))>tol
        error("Test_polyopvar_vector_states failed: scalar combination mismatch (trial=%d).",trial)
    end

    % '+' between two INDEPENDENTLY-valued, identically-shaped state vectors
    % (same dims_in/vars/dom, different names and values)
    y = x;
    y.varname = cellfun(@(s)['y',s(2:end)],x.varname,'UniformOutput',false);
    yvals = cell(p,1);
    for i=1:p
        yvals{i} = rand_poly([dims_in(i),1],vars(depmat_in(i,:)),2,5);
    end
    z4 = x+y;
    Y4 = apply_poly(z4,[xvals;yvals]);
    err4 = cleanpoly(Y4-(expected+vertcat(yvals{:})),1e-10);
    if max(max(abs(err4.C)))>tol
        error("Test_polyopvar_vector_states failed: x+y mismatch for vector states (trial=%d).",trial)
    end
end
disp('Test_polyopvar_vector_states Part 1 (arith) passed.')

%% Part 2: combine_terms on a vector-valued (dims.in>1) single state
for trial = 1:20
    rng(12000+trial);

    [w,~,wdims,wvars] = rand_polyopvar(1,3);  % single state, size in {1,2,3}
    w2 = w.*w;
    if size(w2.Z,1)~=1
        error("Test_polyopvar_vector_states failed: setup assumption violated (trial=%d).",trial)
    end

    % Duplicate the single monomial row (and its coefficient block) by hand,
    % exactly as Test_polyopvar_combine_terms.m does for the scalar case
    w3 = w2;
    w3.Z = [w2.Z; w2.Z];
    w3.C.ops = [w2.C.ops, w2.C.ops];
    w3.C.depmat.in = [w2.C.depmat.in; w2.C.depmat.in];
    w3.C.dims.in = [w2.C.dims.in; w2.C.dims.in];
    w3.chkval = false; w3.C.chkval = false;

    w4 = combine_terms(w3);
    if size(w4.Z,1)~=1
        error("Test_polyopvar_vector_states failed: duplicate rows not merged for vector state (trial=%d,d=%d).", ...
            trial,wdims(1))
    end

    wval = rand_poly([wdims(1),1],wvars,2,4);
    Yorig = apply_poly(w2,{wval});
    Ydup  = apply_poly(w3,{wval});
    Ymerged = apply_poly(w4,{wval});
    errdup = cleanpoly(Ydup-2*Yorig,1e-10);
    if max(max(abs(errdup.C)))>tol
        error("Test_polyopvar_vector_states failed: manual duplication setup did not read as 2x (trial=%d).",trial)
    end
    errmerged = cleanpoly(Ymerged-Ydup,1e-10);
    if max(max(abs(errmerged.C)))>tol
        error("Test_polyopvar_vector_states failed: combine_terms changed the value for a vector state (trial=%d).",trial)
    end
end
disp('Test_polyopvar_vector_states Part 2 (combine_terms) passed.')

%% Part 3: 3-5 states, MIXED (scalar and vector) sizes, vertcat/size
for trial = 1:20
    rng(13000+trial);

    p = randi([3,5]);
    dims_in = randi([1,3],p,1);
    [x,~,dims_in,vars,dom,depmat_in] = rand_polyopvar(p,dims_in);

    xvals = cell(p,1);
    for i=1:p
        xvals{i} = rand_poly([dims_in(i),1],vars(depmat_in(i,:)),2,4);
    end
    Y = apply_poly(x,xvals);
    expected = vertcat(xvals{:});
    err = cleanpoly(Y-expected,1e-10);
    if max(max(abs(err.C)))>tol
        error("Test_polyopvar_vector_states failed: %d-state mixed-size vertcat mismatch (trial=%d,dims_in=%s).", ...
            p,trial,mat2str(dims_in'))
    end
    if ~isequal(size(x),[sum(dims_in),1])
        error("Test_polyopvar_vector_states failed: %d-state mixed-size size() mismatch (trial=%d).",p,trial)
    end
end
disp('Test_polyopvar_vector_states Part 3 (3+ mixed-size states) passed.')

disp('Test_polyopvar_vector_states passed.')
