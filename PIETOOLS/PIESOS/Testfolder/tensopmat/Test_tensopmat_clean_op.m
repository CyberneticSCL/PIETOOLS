% Test_tensopmat_clean_op.m -- tests that @tensopmat/clean_op.m cleans
% every (nonempty) block by dispatching to @tensopvar/clean_op.m,
% correctly determines whether the WHOLE block matrix is zero, and never
% empties a zeroed block (since Cop.ops{i,j}=[] is reserved for "block not
% provided", and initialize.m requires at least one nonempty block per row
% and per column). Swept over varying numbers of row/column blocks and
% spatial variables, and over single- vs multi-term blocks.
% Oracle: apply_op of the cleaned vs. the uncleaned operator on the same
% random input.

% DJ (CC), 10/15/2026: Initial coding;

tol = 1e-8;

for trial = 1:10
    rng(6600+trial);

    p = randi(3)+1;   d = randi(3)+1;
    N = randi(3)+1;   M = randi(3)+1;
    nterms_max = 1 + mod(trial,2);   % alternate single-term / multi-term blocks

    Pop = rand_tensopmat(p,d,N,M,[0,1],0.7,0.8,struct(),'s','t',1,nterms_max,2);
    [p,d] = size(Pop.ops);
    populated = find(~cellfun(@isempty,Pop.ops));
    if numel(populated)<2
        continue   % need at least 2 populated blocks to test a non-trivial drop
    end

    % Input covering every column, from the first nonempty block in each
    X = cell(1,d);
    for j=1:d
        i0 = find(~cellfun(@isempty,Pop.ops(:,j)),1);
        vars_ij = Pop.ops{i0,j}.vars.in;
        X{j} = rand_poly([Pop.dims.in(j),1],polynomial(vars_ij(:)),2,20);
    end

    %%%% Case 1: zero out exactly one populated block -- not flagged zero,
    %%%% action unaffected.
    idx1 = populated(1);
    [i1,j1] = ind2sub([p,d],idx1);
    Pop_bad = Pop;
    Pop_bad.ops{i1,j1} = 0*Pop_bad.ops{i1,j1};
    [Pop_clean,is_zero1] = clean_op(Pop_bad);
    if is_zero1
        error("Test_tensopmat_clean_op failed: operator incorrectly flagged entirely zero after zeroing only one block (trial=%d).",trial)
    elseif isempty(Pop_clean.ops{i1,j1})
        error("Test_tensopmat_clean_op failed: a zeroed block was emptied out (trial=%d).",trial)
    end
    Y_bad = apply_op(Pop_bad,X);
    Y_clean = apply_op(Pop_clean,X);
    err = cleanpoly(Y_bad-Y_clean,1e-10);
    if max(max(abs(full(err.C))))>tol
        error("Test_tensopmat_clean_op failed: cleaned operator does not act identically (trial=%d).",trial)
    end

    %%%% Case 2: zero out EVERY populated block -- is_zero true, zero map.
    Pop_allzero = Pop;
    for k=1:numel(populated)
        [ik,jk] = ind2sub([p,d],populated(k));
        Pop_allzero.ops{ik,jk} = 0*Pop_allzero.ops{ik,jk};
    end
    [Pop_allzero_clean,is_zero2] = clean_op(Pop_allzero);
    if ~is_zero2
        error("Test_tensopmat_clean_op failed: entirely-zero operator not flagged is_zero (trial=%d).",trial)
    end
    for k=1:numel(populated)
        [ik,jk] = ind2sub([p,d],populated(k));
        if isempty(Pop_allzero_clean.ops{ik,jk})
            error("Test_tensopmat_clean_op failed: a zeroed block was emptied out in the all-zero case (trial=%d).",trial)
        end
    end
    Y_allzero = apply_op(Pop_allzero_clean,X);
    if max(max(abs(full(double(cleanpoly(Y_allzero,1e-10).C)))))>tol
        error("Test_tensopmat_clean_op failed: all-zero case did not act as the zero operator (trial=%d).",trial)
    end

    %%%% Case 3: nothing zero -- is_zero false, action unchanged.
    [Pop_clean3,is_zero3] = clean_op(Pop);
    if is_zero3
        error("Test_tensopmat_clean_op failed: nonzero operator incorrectly flagged is_zero (trial=%d).",trial)
    end
    Y0 = apply_op(Pop,X);
    Y3 = apply_op(Pop_clean3,X);
    err3 = cleanpoly(Y0-Y3,1e-10);
    if max(max(abs(full(err3.C))))>tol
        error("Test_tensopmat_clean_op failed: cleaning an already-clean operator changed its action (trial=%d).",trial)
    end

    %%%% Case 4: sub-tolerance noise in one block's sopvar coefficient is
    %%%% scrubbed, other blocks untouched.
    ztol = 1e-6;
    noise = ztol/1000;
    idx4 = populated(end);
    [i4,j4] = ind2sub([p,d],idx4);
    Pop_noisy = Pop;
    Pop_noisy.ops{i4,j4}.ops{1,1}.params{1}(1,1) = noise;
    [Pop_clean4,~] = clean_op(Pop_noisy,ztol);
    if Pop_clean4.ops{i4,j4}.ops{1,1}.params{1}(1,1)~=0
        error("Test_tensopmat_clean_op failed: sub-tolerance noise in a block was not zeroed (trial=%d).",trial)
    end
    for k=1:numel(populated)
        if populated(k)==idx4
            continue
        end
        [ik,jk] = ind2sub([p,d],populated(k));
        if ~isequal(Pop_clean4.ops{ik,jk},Pop.ops{ik,jk})
            error("Test_tensopmat_clean_op failed: cleaning noise in one block changed an untouched block (trial=%d).",trial)
        end
    end
end

disp('Test_tensopmat_clean_op passed.')
